function [fitness,output_vector,d_conv] = mandle_fitness(res, dims_x, input_vector,image_target)

extent = 9;

a_r_s = interp1([0,1],[0,10],input_vector(1));
a_r_e = interp1([0,1],[-extent,extent],input_vector(2));
a_i_s = interp1([0,1],[0,10],input_vector(3));
a_i_e = interp1([0,1],[-extent,extent],input_vector(4));

c_r_s = interp1([0,1],[0,10],input_vector(5));
c_r_e = interp1([0,1],[-extent,extent],input_vector(6));
c_i_s = interp1([0,1],[0,10],input_vector(7));
c_i_e = interp1([0,1],[-extent,extent],input_vector(8));

t1_r_s = interp1([0,1],[0,10],input_vector(9));
t1_r_e = interp1([0,1],[-extent,0],input_vector(10));
t1_i_s = interp1([0,1],[0,10],input_vector(11));
t1_i_e = interp1([0,1],[-extent,0],input_vector(12));

t2_r_s = interp1([0,1],[0,10],input_vector(13));
t2_r_e = interp1([0,1],[-extent,extent],input_vector(14));
t2_i_s = interp1([0,1],[0,10],input_vector(15));
t2_i_e = interp1([0,1],[-extent,extent],input_vector(16));

a = a_r_s*10^a_r_e + j*a_i_s*10^a_i_e;
c = c_r_s*10^c_r_e + j*c_i_s*10^c_i_e;
t1 = t1_r_s*10^t1_r_e + j*t1_i_s*10^t1_i_e;
t2 = t2_r_s*10^t2_r_e + j*t2_i_s*10^t2_i_e;

output_vector=[
a
c
t1
t2
];

d_conv = run_mandle(res, dims_x, a, c, t1, t2);

conv_resized = imresize(d_conv,size(image_target),"nearest");
conv_region = conv_resized==0;

fitness = 1/immse(uint8(image_target), uint8(conv_region));

fitness = fitness + log10(numel(unique(d_conv)));

% subplot(1,2,1)
% imagesc(d_conv)
% axis tight equal
% 
% subplot(1,2,2)
% imagesc(conv_region)
% axis tight equal

end
