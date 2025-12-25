

%%%%%%% Code to generate Sparse Quadrature Points %%%%%%%%%

uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);


[status,out] = system(['./generate_quad -d ', num2str(dim), ' -p ', num2str(num_qd), ' -g HG -x full'])

% load the quadrature points and weights
xi = load('qdpts.dat'); % Need to scale each variable, qdpts = f(qdpts)
w = load('wghts.dat');
 