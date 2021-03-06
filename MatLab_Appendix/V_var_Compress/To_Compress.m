
I = imread('barnsley.jpg');
I = imresize(I, [512,512]);
I = rgb2gray(I);
figure(1)
imshow(I)

i = 0;

[j,k] = Fractal_VC(I,i);


%M = [[1,1,4,3,3,3,2,1];[3,2,1,2,3,3,4,1];[1,1,4,3,3,3,2,1];[3,2,1,2,3,3,4,1];[2,2,2,4,1,4,3,2];[2,2,3,4,1,4,3,2];[2,2,2,4,1,4,3,2];[4,2,3,4,1,4,3,1];[1,3,1,2,1,1,1,3];[3,4,1,2,1,1,1,3];[3,3,1,2,1,1,1,3];[2,4,1,2,1,1,1,3];[2,4,3,1,2,2,3,4];[4,2,2,1,4,2,3,4];[4,3,3,1,2,2,3,4];[4,3,2,1,4,2,3,4]];
%MQ = [[138];[33];[171];[37]];

reco = Decode(j,k);
figure(2)
fig = imshow(reco,[0,256])
saveas(fig,string(i)+'.png');
