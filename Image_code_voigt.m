f = @(x, y, z0, A, mu, w1, w2, xc, yc) ...
    z0 + A * ( ...
    ((mu * 4 / pi^2) * w1 * w2) ./ (1 + 4 * ((x - xc) ./ w1).^2 + 4 * ((y - yc) ./ w2).^2) ...
    + (1 - mu) * (4 * log(2)) / (pi * w1 * w2) * exp(-4 * log(2) / w1^2 * (x - xc).^2 - 4 * log(2) / w2^2 * (y - yc).^2) );

%electron repulsion term
vx_rand=rep_std*randn(MAX,1); %sigma = rep_std를 가지는 normaldistribution 해당하는 MAX개 추출
vy_rand=rep_std*randn(MAX,1);
vz_rand=rep_std*randn(MAX,1);
Cell_velocity_el=cell(numel(I2_list),1);
vx_pix=cell(numel(I2_list),1);
vy_pix=cell(numel(I2_list),1);

v_rand=sqrt(vx_rand.^2 + vy_rand.^2 + vz_rand.^2 ); %위에서 뽑은 입자들의 속도 크기계산
for ind = 1:numel(I2_list)
Cell_velocity_el{ind}(:,1)=Cell_velocity{ind}(:,1)+rep_vel*vx_rand./v_rand; %x방향으로 project시킨 뒤 rep_vel을 곱해서 recoil velocity줌.
Cell_velocity_el{ind}(:,2)=Cell_velocity{ind}(:,2)+rep_vel*vy_rand./v_rand;
vx_pix{ind} = Cell_velocity_el{ind}(:,1)./vel_pix;
vy_pix{ind} = Cell_velocity_el{ind}(:,2)./vel_pix;
end

%image blur term
xgrid=-125:125;
ygrid=-125:125;

[xgird_2, ygrid_2]=meshgrid(xgrid,ygrid);
z0=-1.5; A=1117567;
xc=0; yc=0; mu=6.3e-4;
w1=6.306*0.92; w2=8.015*0.92;
M=f(xgird_2, ygrid_2, z0, A, mu, w1, w2, xc, yc);

h2d=cell(numel(I2_list),1);
img_matrix=cell(numel(I2_list),1);
for ind = 1:numel(I2_list)
    h2d{ind}=hist3([vx_pix{ind}(:,1),vy_pix{ind}(:,1)],'Ctrs',{xgrid ygrid});
    img_matrix{ind}=conv2(h2d{ind},M,'same');
end

sumx=zeros(numel(xgrid),numel(I2_list));
sumy=zeros(numel(ygrid),numel(I2_list));

for ind=1:numel(I2_list)
sumx(:,ind)=sum(img_matrix{ind},2);
sumy(:,ind)=sum(img_matrix{ind},1);
end