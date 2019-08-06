%%%% demo code  'Haze removal for a single visible remote sensing image'
%%%% 20131206.tif: test landsat8 image
%%%% guidedfilter.m,boxfilter.m: guidedfilter
%%%% rsdehaze.p: main function
%%%% rsshow.m: show rs image using linear stretch
%%%% TVinpaint.m,tv_inpaint.mexw64: we use a fast image inpainting tool in  http://www.imm.dtu.dk/~pcha/mxTV/
clear,clc;

ms = imread( '20131206.tif' );
%%% landsat8 B,G,R band, select correct band for other satellite images
for i=1:3
    band{i}=double( ms(:,:,i+1) );
end
figure, rsshow(band{1},band{2},band{3} ,1  );

bande = band{1} + band{1}- 0.9*band{2};% Create synthetic band
bande = max(bande,0);
maskth = 0.001; % threshold for brightest pixels
par_win = 10; % patch rad 
par_gf1 = 0.1; % GF filtering eps
par_gf2 = 2*par_win-1; % GF filtering rad
bande_edge = edge(bande,'canny',0.1);
out = rsdehaze(band,bande_edge,maskth, par_win,par_gf1,par_gf2);
figure,rsshow(out{1},out{2},out{3},1);


