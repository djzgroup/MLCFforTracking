
%PRECISION_PLOT
%   Calculates precision for a series of distance thresholds (percentage of
%   frames where the distance to the ground truth is within the threshold).
%   The results are shown in a new figure if SHOW is true.
%
%   Accepts positions and ground truth as Nx2 matrices (for N frames), and
%   a title string.
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/

%%%%%%%%%%%%%%% plotting  precision picture of SRDCF%%%%%%%%%%%%%%%%%%%%%
fid_1 = fopen('test.txt','r');
fid_2 = fopen('truth.txt','r');
A = zeros(140,4);
B = zeros(140,4);

for i=1:140
    for j=1:3
        A(i,j) = fscanf(fid_1,'%f',[1,1]);B(i,j) = fscanf(fid_2,'%f',[1,1]);
    end
        A(i,4) = fscanf(fid_1,'%f',[1,1]);B(i,4) = fscanf(fid_2,'%f',[1,1]);
end
fclose(fid_1);
fclose(fid_2);
clear fid_1 fid_2;
%precisions = precision_plot(A, B, 'sequences/Couple', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_3 = fopen('test.txt','r');
fid_4 = fopen('truth.txt','r');
C = zeros(140,4);
D = zeros(140,4);

for i=1:140
    for j=1:3
        C(i,j) = fscanf(fid_3,'%f',[1,1]);D(i,j) = fscanf(fid_4,'%f',[1,1]);
    end
        C(i,4) = fscanf(fid_3,'%f',[1,1]);D(i,4) = fscanf(fid_4,'%f',[1,1]);
end
fclose(fid_3);
fclose(fid_4);
clear fid_3 fid_4;
%precisions = precision_plot(A, B, 'sequences/Couple', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	max_threshold = 50;  %used for graphs in the paper
	positions_1 = A;
	positions_2 = C;
	ground_truth = B;
	
	precisions_1 = zeros(max_threshold, 1);
	precisions_2 = zeros(max_threshold, 1);
	
	if size(positions_1,1) ~= size(ground_truth,1),
% 		fprintf('%12s - Number of ground truth frames does not match number of tracked frames.\n', title)
		
		%just ignore any extra frames, in either results or ground truth
		n = min(size(positions_1,1), size(ground_truth,1));
		positions_1(n+1:end,:) = [];
		positions_2(n+1:end,:) = [];
		ground_truth(n+1:end,:) = [];
	end
	
	%calculate distances to ground truth over all frames
	distances_1 = sqrt((positions_1(:,1) - ground_truth(:,1)).^2 + ...
				 	 (positions_1(:,2) - ground_truth(:,2)).^2);
	distances_1(isnan(distances_1)) = [];
	
	distances_2 = sqrt((positions_2(:,1) - ground_truth(:,1)).^2 + ...
				 	 (positions_2(:,2) - ground_truth(:,2)).^2);
	distances_2(isnan(distances_2)) = [];
	
    aaa = 0;
	for kk = 1:size(distances_1,1)
		aaa = aaa + distances_2(kk);
	end
	aaa = aaa/kk;
	clear kk;
	disp('there CLE :');
	disp(aaa);
    
    
	precisions_U = A(:,1:2) + A(:,3:4);
	ground_truth_U = B(:,1:2) + B(:,3:4);
	%A B
	precision_overlap_min_1 = min([precisions_U(:,1),ground_truth_U(:,1)],[],2);
	precision_overlap_min_2 = min([precisions_U(:,2),ground_truth_U(:,2)],[],2);
	precision_overlap_min = [precision_overlap_min_1,precision_overlap_min_2];
	
	precision_overlap_max_1 = max([A(:,1),B(:,1)],[],2);
	precision_overlap_max_2 = max([A(:,2),B(:,2)],[],2);
	precision_overlap_max = [precision_overlap_max_1,precision_overlap_max_2];
	
	clear precision_overlap_max_1 precision_overlap_max_2 precision_overlap_min_1 precision_overlap_min_2;
	
	size_overlap = abs(precision_overlap_min - precision_overlap_max);
	%size_no_overlap = abs(precision_no_overlap_max - precision_no_overlap_min);
	
	Area_overlap = size_overlap(:,1) .* size_overlap(:,2);
	Area_no_overlap = ((A(:,3) .* A(:,4)) + (B(:,3) .* B(:,4))) - Area_overlap;
	
	precision_overlap = Area_overlap ./ Area_no_overlap;
	precis_overlap = zeros(100,1);
	for pp = 1:100,
		precis_overlap(pp) = nnz(precision_overlap <= (pp/100)) / numel(precision_overlap);
	end
	
	figure('Number','off', 'Name',['Precision_overlap - ' ' precision_overlap'])
	plot(precis_overlap, 'k-', 'LineWidth',2)
	xlabel('Threshold'), ylabel('Precision_overlap')
	
	%compute precisions
	for p = 1:max_threshold,
		precisions_2(p) = nnz(distances_2 <= p) / numel(distances_2);
    end
    
	for p = 1:max_threshold,
		precisions_1(p) = nnz(distances_1 <= p) / numel(distances_1);
	end
	%plot the precisions
	
		figure('Number','off', 'Name',['Precisions - ' 'one_1'])
		plot(precisions_1, 'r', 'LineWidth',2)
        hold on;
		plot(precisions_2, 'k-', 'LineWidth',2)
		xlabel('Threshold'), ylabel('Precision')

	


