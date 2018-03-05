function saveQuad(working_folder,quadData,allThresholds,allHeights,uvHisData,edgesAll,pdfHeights)
save([ working_folder '/quadData'],'quadData','allThresholds','allHeights')
save([ working_folder '/histData'],'uvHisData','edgesAll','pdfHeights')
end

