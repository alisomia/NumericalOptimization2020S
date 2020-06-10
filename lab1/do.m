ind = 1;
% prob = {die2 die10 die20 die30 die40 die50};
% init = {dieinit(2) dieinit(10) dieinit(20) dieinit(30) dieinit(40) dieinit(50)};
prob = {mse3 mse5 mse7};
init = {mseinit(3) mseinit(5) mseinit(7)};
% prob = {bd4, bd10, bd20, bd30, bd40, bd50};
% init = {[25;5;-5;1] [25;5;-5;1] [25;5;-5;1] [25;5;-5;1] [25;5;-5;1] [25;5;-5;1]};
for i = 1:3 %[bd4 bd10, bd20,bd30,bd40,bd50]
[x, info] = quasiNewton(prob{i}, init{i},[],"quad",[],10000,[],"mixed",0);
  RESULT(ind, 1:5) = [info.iter, info.subiter, info.time, info.call, info.f];
  [x, info] = quasiNewton(prob{i}, init{i},[],"quad",[],10000,[],"mixed",1/2);
    RESULT(ind,6:10) = [info.iter, info.subiter, info.time, info.call, info.f];
  [x, info] = quasiNewton(prob{i}, init{i},[],"quad",[],10000,[],"mixed",1);
    RESULT(ind + 4, 1:5) = [info.iter, info.subiter, info.time, info.call, info.f];
  [x, info] = quasiNewton(prob{i}, init{i},[],"quad",[],10000,[],"mixed","SR1");
    RESULT(ind + 4, 6:10) = [info.iter, info.subiter, info.time, info.call, info.f];
    ind = ind + 1;
end