// pybind11 binding for the vendored realSEUDO Fista::Seudo solver
// (see vendor/LICENSE). This is a fresh binding (not derived from the
// MATLAB MEX gateway vendor/../seudo_native_newapi.cpp, though it mirrors
// its parameter interface/semantics), targeting Python via numpy arrays
// instead of MATLAB mxArrays.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>

#include "vendor/fista.hpp"
#include "vendor/fista_gradient.hpp"
#include "vendor/fista_types.hpp"
#include "vendor/seudo.hpp"

namespace py = pybind11;

// Inputs:
//   image[ht][wd]        - frame/window to fit, row-major (numpy default)
//   blob[blobht][blobwd] - blob basis function, row-major
//   blob_spacing         - >=1, pixel spacing of the blob grid (1 = one blob per pixel)
//   rois[wd*ht][n_rois]  - each column is one profile, flattened row-major
//                          to match `image`'s (y*wd+x) pixel ordering
//   weights0             - initial weights (n_rois + wd*ht), or empty for all-0s
//   lambda_               - per-weight L1 lambda (n_rois + wd*ht), or empty for all-0s
//   eps, max_steps        - FISTA stopping tolerance / iteration cap
//   l_mode                - 0: dynamic L, 1: static multi-L, 2: dynamic L + fast brake,
//                           3: static multi-L + fast brake (2 matches this project's Python solver)
//   stop_mode              - 0: relative norm2 (matches TFOCS / the Python solver), 1: norm2,
//                            2: every dimension
//   verbose, nthreads
//
// Returns (weights, n_steps, log).
static py::tuple seudo_solve(
    py::array_t<double, py::array::c_style | py::array::forcecast> image,
    py::array_t<double, py::array::c_style | py::array::forcecast> blob,
    double blob_spacing,
    py::array_t<double, py::array::c_style | py::array::forcecast> rois,
    py::array_t<double, py::array::c_style | py::array::forcecast> weights0,
    py::array_t<double, py::array::c_style | py::array::forcecast> lambda_in,
    double eps,
    int max_steps,
    int l_mode,
    int stop_mode,
    bool verbose,
    int nthreads
) {
    auto image_buf = image.request();
    if (image_buf.ndim != 2)
        throw std::invalid_argument("image must be a 2D (ht, wd) array");
    int ht = static_cast<int>(image_buf.shape[0]);
    int wd = static_cast<int>(image_buf.shape[1]);

    auto seudo = std::make_shared<Fista::Seudo>(wd, ht);
    {
        const double *p = static_cast<const double *>(image_buf.ptr);
        seudo->image_.img_.assign(p, p + static_cast<size_t>(wd) * ht);
    }

    auto blob_buf = blob.request();
    if (blob_buf.ndim != 2)
        throw std::invalid_argument("blob must be a 2D array");
    int blobht = static_cast<int>(blob_buf.shape[0]);
    int blobwd = static_cast<int>(blob_buf.shape[1]);
    auto blob_sprite = std::make_shared<Fista::Sprite>(blobwd, blobht, blobwd / 2, blobht / 2);
    {
        const double *p = static_cast<const double *>(blob_buf.ptr);
        blob_sprite->img_.assign(p, p + static_cast<size_t>(blobwd) * blobht);
    }
    seudo->blob_ = blob_sprite;

    if (blob_spacing < 1.0)
        throw std::invalid_argument("blob_spacing must be >= 1");
    seudo->blobSpacing_ = static_cast<int>(blob_spacing);

    auto rois_buf = rois.request();
    if (rois_buf.ndim != 2)
        throw std::invalid_argument("rois must be a 2D (wd*ht, n_rois) array");
    int rois_d0 = static_cast<int>(rois_buf.shape[0]);
    int n_rois = static_cast<int>(rois_buf.shape[1]);
    if (rois_d0 != wd * ht)
        throw std::invalid_argument("rois first dimension must equal wd*ht of image");
    {
        const double *p = static_cast<const double *>(rois_buf.ptr);
        for (int r = 0; r < n_rois; r++) {
            auto sprite = std::make_shared<Fista::Sprite>(wd, ht, 0, 0);
            for (int pix = 0; pix < wd * ht; pix++) {
                sprite->img_[pix] = p[static_cast<size_t>(pix) * n_rois + r];
            }
            seudo->rois_.emplace_back(sprite);
        }
    }

    int n_weights_total = n_rois + wd * ht;

    auto weights_buf = weights0.request();
    size_t n_w = static_cast<size_t>(weights_buf.size);
    if (n_w == 0) {
        seudo->weights_.assign(n_weights_total, 0.0);
    } else {
        if (static_cast<int>(n_w) != n_weights_total)
            throw std::invalid_argument("weights0 must have (n_rois + wd*ht) elements, or be empty");
        const double *p = static_cast<const double *>(weights_buf.ptr);
        seudo->weights_.assign(p, p + n_w);
    }

    auto lambda_buf = lambda_in.request();
    size_t n_l = static_cast<size_t>(lambda_buf.size);
    if (n_l == 0) {
        seudo->lambda_.assign(n_weights_total, 0.0);
    } else {
        if (static_cast<int>(n_l) != n_weights_total)
            throw std::invalid_argument("lambda_ must have (n_rois + wd*ht) elements, or be empty");
        const double *p = static_cast<const double *>(lambda_buf.ptr);
        seudo->lambda_.assign(p, p + n_l);
    }

    seudo->eps_ = eps;
    seudo->maxSteps_ = max_steps;

    switch (l_mode) {
        case 0: seudo->multiGrad_ = false; seudo->fastBrake_ = false; break;
        case 1: seudo->multiGrad_ = true;  seudo->fastBrake_ = false; break;
        case 2: seudo->multiGrad_ = false; seudo->fastBrake_ = true;  break;
        case 3: seudo->multiGrad_ = true;  seudo->fastBrake_ = true;  break;
        default: throw std::invalid_argument("l_mode must be 0-3");
    }

    switch (stop_mode) {
        case 0: seudo->stopping_ = Fista::Run::StopEpsNorm2Rel; break;
        case 1: seudo->stopping_ = Fista::Run::StopEpsNorm2; break;
        case 2: seudo->stopping_ = Fista::Run::StopEpsEveryDimension; break;
        default: throw std::invalid_argument("stop_mode must be 0-2");
    }

    seudo->verbose_ = verbose;
    seudo->setNumThreads(nthreads);

    {
        // release the GIL for the duration of the (potentially multithreaded) solve
        py::gil_scoped_release release;
        seudo->compute();
    }

    if (!seudo->error_.empty())
        throw std::runtime_error(seudo->error_);

    py::array_t<double> out_weights(n_weights_total);
    {
        auto out_buf = out_weights.request();
        double *op = static_cast<double *>(out_buf.ptr);
        for (int i = 0; i < n_weights_total; i++) op[i] = seudo->weights_[i];
    }

    return py::make_tuple(out_weights, seudo->stepsTaken_, seudo->log_);
}

PYBIND11_MODULE(_seudo_native, m) {
    m.doc() = "Native (C++) accelerator for the SEUDO per-frame FISTA solve, "
              "wrapping realSEUDO's Fista::Seudo class (vendored, MIT licensed, "
              "see vendor/LICENSE).";
    m.def("seudo_solve", &seudo_solve,
        py::arg("image"), py::arg("blob"), py::arg("blob_spacing"),
        py::arg("rois"), py::arg("weights0"), py::arg("lambda_"),
        py::arg("eps"), py::arg("max_steps"), py::arg("l_mode"),
        py::arg("stop_mode"), py::arg("verbose"), py::arg("nthreads"),
        "Run the SEUDO FISTA solve for one frame/window; returns (weights, n_steps, log).");
}
