function se = provideOtherClassification(se, newClass)

for cc = 1:se.nCells
    se.tcDefault.transientInfo(cellID).classification            = newClass{cc}.binClass;
    se.tcDefault.transientInfo(cellID).classMetaInfo.profLevel   = newClass{cc}.profLevel;
    se.tcDefault.transientInfo(cellID).classMetaInfo.blobLevels  = newClass{cc}.blobLevels;
end

end
