rule.name ="wolfe";
rule.rho = 0.3;
rule.sigma = 0.6;
%[~,info] = conjGrad(p, [], [], [], method, 3000, 0, 1e-6, 0)
%[~,info] = conjGrad(ep4,[],[],[],"FR",20000,5,1e-6,20)
%[~,info] = conjGrad(trid4,[],[],[],"PRP+",20000,5,1e-6,0)
%[~,info] = conjGrad(mat4,[],[],[],"FR-PRP",60000,5,1e-6,50)
%[~,info] = BB(mat3,[], [] ,1000000, 1e-6)
RESULT = zeros(12,5);
indx = 1;
global cf cg
cf = 0;
cg = 0;
[~,info] = conjGrad(mat2,[],[],[],100,100000,5,1e-6,0)

% for p = [mat2, mat3]%[trid2, trid3, trid4] %[ep2, ep3, ep4]%trig2, trig3, trig4]
% 
%     indy = 1;
%     for method = ["FR","PRP+","FR-PRP",]
%         cf = 0; cg = 0;
%         [~,info] = conjGrad(p,[],[],[],method,60000,5,1e-6,50)
%         RESULT(indx,indy) = info.time;
%         RESULT(indx+1,indy) = cf;
%         RESULT(indx+2, indy) = cg;
%         RESULT(indx+3, indy) = info.iter;
%         indy = indy+1;
%     end
% %     cf = 0; cg = 0;
% %     [~,info] = BB(p,[], [] ,10000, 1e-6)
% %         RESULT(indx,indy) = info.time;
% %         RESULT(indx+1,indy) = cf;
% %         RESULT(indx+2, indy) = cg;
% %         RESULT(indx+3, indy) = info.iter;
%     indx = indx + 4;
% end