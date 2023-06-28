function SlicerPNGs(filenameA, filenameB, prefixA, prefixB, outputdir)
% ClinicalASL toolbox 2023, JCWSiero
    %create Slicer PNGs for outlines of two overlapping dataset
    currentdir = pwd;
    cd(outputdir)
    eval(['!slicer ' filenameA ' ' filenameB ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -y 0.75 slm.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']);
    eval(['!pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + slm.png + sli.png + slj.png + slk.png + sll.png ' outputdir '/' prefixA '2' prefixB '1.png ']);
    eval(['!slicer ' filenameB ' ' filenameA ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -y 0.75 slm.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']);
    eval(['!pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + slm.png + sli.png + slj.png + slk.png + sll.png ' outputdir '/' prefixA '2' prefixB '2.png ']);
    eval(['!pngappend ' outputdir '/' prefixA '2' prefixB '1.png - ' outputdir '/' prefixA '2' prefixB '2.png ' outputdir '/' prefixA '2' prefixB '.png']);
    eval(['!/bin/rm -f sl?.png ' outputdir '/' prefixA '2' prefixB '1.png ' outputdir '/' prefixA '2' prefixB '2.png']);
    cd(currentdir)
