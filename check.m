function []=check(aaa,bbb)
						disp(bbb);				%%
	disp('function check parameters');
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				aa=input('inpust a number:','s');%%
				if aa=='1'						%%
					disp(size(aaa));			%%
				elseif aa=='2'					%%
					disp(aaa);					%%
					disp(size(aaa));			%%
				elseif aa=='2'					%%
					h=figure(852);				%%
					xxx=[1:size(aaa,1)];		%%
					yyy=[1:size(aaa,2)];		%%
					mesh(xxx,yyy,aaa);			%%	
					surf(xxx,yyy,aaa);			%%
					disp('input any key to leave!!');
					pause;
					close(h);					%%	
				else							%%
					h=figure(852);				%%
					xxx=[1:size(aaa,1)];		%%
					yyy=[1:size(aaa,2)];		%%
					mesh(xxx,yyy,aaa);			%%	
					surf(xxx,yyy,aaa);			%%
					disp(aaa);					%%
					disp('input any key to leave!!');
					pause;						%%
					close(h);					%%	
				end								%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('leave check function!');
end