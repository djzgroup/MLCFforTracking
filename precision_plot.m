function precisions = precision_plot(positions, ground_truth, title, show)
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

	
	max_threshold = 50;  %used for graphs in the paper
	
	
	precisions = zeros(max_threshold, 1);
	
	if size(positions,1) ~= size(ground_truth,1),
% 		fprintf('%12s - Number of ground truth frames does not match number of tracked frames.\n', title)
		
		%just ignore any extra frames, in either results or ground truth
		n = min(size(positions,1), size(ground_truth,1));
		positions(n+1:end,:) = [];
		ground_truth(n+1:end,:) = [];
	end
	
	%calculate distances to ground truth over all frames
	distances = sqrt((positions(:,1) - ground_truth(:,1)).^2 + ...
				 	 (positions(:,2) - ground_truth(:,2)).^2);
	distances(isnan(distances)) = [];

	%%%%%%%%%%%%	CLE		%%%%%%%%%%%%%%
	aaa = 0;
	for kk = 1:size(distances,1)
		aaa = aaa + distances(kk);
	end
	aaa = aaa/kk;
	clear kk;
	disp('there CLE :');
	disp(aaa);
	%%%%%%%%%%%%%	OP		%%%%%%%%%%%%%%%%%%%%%%%%
	precisions_U = positions(:,1:2) + positions(:,3:4);
	ground_truth_U = ground_truth(:,1:2) + ground_truth(:,3:4);
	
		precision_overlap_min_1 = min([precisions_U(:,1),ground_truth_U(:,1)],[],2);
	precision_overlap_min_2 = min([precisions_U(:,2),ground_truth_U(:,2)],[],2);
	precision_overlap_min = [precision_overlap_min_1,precision_overlap_min_2];
	
	precision_overlap_max_1 = max([positions(:,1),ground_truth(:,1)],[],2);
	precision_overlap_max_2 = max([positions(:,2),ground_truth(:,2)],[],2);
	precision_overlap_max = [precision_overlap_max_1,precision_overlap_max_2];
	
	clear precision_overlap_max_1 precision_overlap_max_2 precision_overlap_min_1 precision_overlap_min_2;
	
	size_overlap = abs(precision_overlap_min - precision_overlap_max);
	%size_no_overlap = abs(precision_no_overlap_max - precision_no_overlap_min);
	
	Area_overlap = size_overlap(:,1) .* size_overlap(:,2);
	Area_no_overlap = ((positions(:,3) .* positions(:,4)) + (ground_truth(:,3) .* ground_truth(:,4))) - Area_overlap;
	
	precision_overlap = Area_overlap ./ Area_no_overlap;
	precis_overlap = zeros(100,1);
	for pp = 1:100,
		precis_overlap(pp) = nnz(precision_overlap <= (pp/100)) / numel(precision_overlap);
	end
	
	
	
	%compute precisions
	for p = 1:max_threshold,
		precisions(p) = nnz(distances <= p) / numel(distances);
	end
	
	%plot the precisions
	if show == 1,
		figure('Number','off', 'Name',['Precisions - ' title])
		plot(precisions, 'k-', 'LineWidth',2)
		xlabel('Threshold'), ylabel('Precision')
		
		figure('Number','off', 'Name',['Precision_overlap - ' title])
		plot(precis_overlap, 'k-', 'LineWidth',2)
		xlabel('Threshold'), ylabel('Precision_overlap')
	end
	
end

