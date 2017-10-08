% for i = 1:Ns
%       ave_err_CMM(i,1)=zeros(1,1);
         var_err_CMM(i)=0;
%         ave_sq_err(i)=0;
%         for k=1:N
%             err_CMM(i,:)=mu'-usrenu(i,1:2);
            ave_err_CMM=err_CMM;
%             bias_x(i)=err_CMM(i,1);
%             bias_y(i)=err_CMM(i,2);
            err_norm(i)=norm(err_CMM(i,:));
            deter(i)=det(cov);
%         end
        for k=1:N
            %var_err_CMM(i)=var_err_CMM(i)+norm((err_CMM{k}(i,:)-ave_err_CMM(i,1:2))')^2/N;
            ave_sq_err(i)=ave_sq_err(i)+norm(err_CMM{k}(i,:)')^2/N;
            norm_ave_err(i)=norm(ave_err_CMM(i,1:2)');
        end
        rt_var_err(i)=sqrt(var_err_CMM(i));
        rt_ave_sq_err(i)=sqrt(ave_sq_err(i));
