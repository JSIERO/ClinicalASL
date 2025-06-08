function SlicerPNGs(filenameA, filenameB, prefixA, prefixB, outputdir)
% ClinicalASL toolbox 2023, JCWSiero
% create Slicer PNGs for outlines of two overlapping dataset
% input filenameA = reference image, filenameB = moved image, prefixA, prefixB= prefixnames for output
% output will be an edge image shwoing overlap 'filenameA' and 'filenameB' with filename 'prefixA2prefixB.png'
currentdir = pwd;
cd(outputdir)
system(['slicer ' filenameA ' ' filenameB ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -y 0.75 slm.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']);
system(['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + slm.png + sli.png + slj.png + slk.png + sll.png ' outputdir '/' prefixA '2' prefixB '1.png ']);
system(['slicer ' filenameB ' ' filenameA ' -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -y 0.75 slm.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png ']);
system(['pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + slm.png + sli.png + slj.png + slk.png + sll.png ' outputdir '/' prefixA '2' prefixB '2.png ']);
system(['pngappend ' outputdir '/' prefixA '2' prefixB '1.png - ' outputdir '/' prefixA '2' prefixB '2.png ' outputdir '/' prefixA '2' prefixB '.png']);
system(['/bin/rm -f sl?.png ' outputdir '/' prefixA '2' prefixB '1.png ' outputdir '/' prefixA '2' prefixB '2.png']);
cd(currentdir)
