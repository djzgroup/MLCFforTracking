function results = OUR_tracker(params)

% parameters
search_area_scale = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
lambda = params.lambda;
learning_rate = params.learning_rate;
refinement_iterations = params.refinement_iterations;
filter_max_area = params.filter_max_area;
nScales = params.number_of_scales;
scale_step = params.scale_step;
interpolate_response = params.interpolate_response;
num_GS_iter = params.num_GS_iter;

%pos_tr = params.pos_tr;

features = params.t_features;

s_frames = params.s_frames;
pos = floor(params.init_pos);
target_sz = floor(params.wsize);

debug = params.debug;
visualization = params.visualization || debug;

num_frames = numel(s_frames);

init_target_sz = target_sz;

%set the feature ratio to the feature-cell size
featureRatio = params.t_global.cell_size;

search_area = prod(init_target_sz / featureRatio * search_area_scale);

% when the number of cells are small, choose a smaller cell size
if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * filter_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));
        
        featureRatio = params.t_global.cell_size;
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area > filter_max_area
    currentScaleFactor = sqrt(search_area / filter_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

%window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;				
rg = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);			
cg = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);			
[rs, cs] = ndgrid( rg,cg);																			
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));												

y1 = circshift(y,[-1,3]);
y2 = circshift(y,[1,-3]);%%%%	1/4 of object_sz 											
y3 = y;																								
yf2 = fft2(y2);																						
yf1 = fft2(y1);																						
yf3 = fft2(y3);																						



if interpolate_response == 1
    interp_sz = use_sz * featureRatio;			
else
    interp_sz = use_sz;
end

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');
%cos_window = ones(size(cos_window));
%cos_window = cos_window*0.8 + 0.2;
%%%there change cos_window
% the search area size
support_sz = prod(use_sz);						

% Calculate feature dimension
im = imread(s_frames{1});
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

% compute feature dimensionality
feature_dim = 0;						
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end;
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end;
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end;
end;

if size(im,3) > 1 && colorImage == false
    im = im(:,:,1);
end

% compute the indices for the real, positive and negative parts of the spectrum

[dft_sym_ind, dft_pos_ind, dft_neg_ind] = partition_spectrum2(use_sz);

% the discrete fourier series output indices
dfs_sym_ind = (1:length(dft_sym_ind))';
dfs_real_ind = dfs_sym_ind(end) - 1 + 2 * (1:length(dft_pos_ind))';
dfs_imag_ind = dfs_sym_ind(end) + 2 * (1:length(dft_pos_ind))';

% construct the transformation matrix from dft to dfs (the real fourier series)
dfs_matrix = dft2dfs_matrix(dft_sym_ind, dft_pos_ind, dft_neg_ind, dfs_sym_ind, dfs_real_ind, dfs_imag_ind);

% create vectorized desired correlation output
yf_vec_1 = single(yf1(:));
yf_vec_2 = single(yf2(:));
yf_vec_3 = single(yf3(:));

if params.use_reg_window
    % create weight window
    ref_window_power = params.reg_window_power;
    
    % normalization factor
    reg_scale = 0.3 * base_target_sz/featureRatio;
    
    % construct grid
    wrg = -(use_sz(1)-1)/2:(use_sz(1)-1)/2;
    wcg = -(use_sz(2)-1)/2:(use_sz(2)-1)/2;
    [wrs, wcs] = ndgrid(wrg, wcg);
    
    % construct the regukarization window
    reg_window = (params.reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^ref_window_power + abs(wcs/reg_scale(2)).^ref_window_power) + params.reg_window_min;

    % compute the DFT and enforce sparsity
    reg_window_dft = fft2(reg_window) / prod(use_sz);
    reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
    reg_window_dft_sep(abs(reg_window_dft_sep) < params.reg_sparsity_threshold * max(abs(reg_window_dft_sep(:)))) = 0;
    reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);
    
    % do the inverse transform, correct window minimum
    reg_window_sparse = real(ifft2(reg_window_dft));
    reg_window_dft(1,1) = reg_window_dft(1,1) - support_sz * min(reg_window_sparse(:)) + params.reg_window_min;
    
    % construct the regularizsation matrix
    regW = cconvmtx2(reg_window_dft);
    
    regW_dfs = real(dfs_matrix * regW * dfs_matrix');
	WW_block = regW_dfs' * regW_dfs;
    % If the filter size is small enough, remove small values in WW_block.
    % It takes too long time otherwise.
    if support_sz <= 120^2
        WW_block(0<abs(WW_block) & abs(WW_block)<0.00001) = 0;
    end
else
    % else use a scaled identity matrix
    WW_block = lambda * speye(support_sz);
    params.reg_window_min = sqrt(lambda);
end

% create block diagonal regularization matrix
WW = eval(['blkdiag(WW_block' repmat(',WW_block', 1, feature_dim-1) ');']);
% upper and lower triangular parts of the regularization matrix
WW_L = tril(WW);
WW_U = triu(WW, 1);

if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    
    scaleFactors = scale_step .^ scale_exp;
    
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

num_sym_coef = length(dft_sym_ind);

% create indexing vectors

% first create the indices for the symmetric (real) part of the spectrum
index_i_sym = zeros(2*feature_dim, length(dft_sym_ind), feature_dim);
index_j_sym = zeros(size(index_i_sym));

index_i_sym_re = repmat(bsxfun(@plus, support_sz*(0:feature_dim-1)', 1:length(dft_sym_ind)), [1 1 feature_dim]); %index for the Real-Real part
index_i_sym(1:2:end, :, :) = index_i_sym_re;
index_i_sym(2:2:end, :, :) = NaN; % these will be zero

index_j_sym_re = permute(index_i_sym_re, [3 2 1]);
index_j_sym(1:2:end, :, :) = index_j_sym_re;
index_j_sym(2:2:end, :, :) = NaN; % these will be zero

% create the indices for the remaining part
index_i = zeros(2*feature_dim, 2*length(dft_pos_ind), feature_dim);
index_j = zeros(size(index_i));


index_i_re = repmat(bsxfun(@plus, support_sz*(0:feature_dim-1)', (length(dft_sym_ind)+1:2:support_sz)), [1 1 feature_dim]); %index for the Real-Real part
index_i(1:2:end, 1:2:end, :) = index_i_re;
index_i(2:2:end, 1:2:end, :) = index_i_re + 1;
index_i(1:2:end, 2:2:end, :) = index_i_re;
index_i(2:2:end, 2:2:end, :) = index_i_re + 1;

index_j_re = permute(index_i_re, [3 2 1]);
index_j(1:2:end, 1:2:end, :) = index_j_re;
index_j(2:2:end, 1:2:end, :) = index_j_re;
index_j(1:2:end, 2:2:end, :) = index_j_re + 1;
index_j(2:2:end, 2:2:end, :) = index_j_re + 1;

% concatenate the results
index_i = cat(2, index_i_sym, index_i);
index_j = cat(2, index_j_sym, index_j);

index_i = index_i(:);
index_j = index_j(:);

% the imaginary part of the autocorrelations (along the diagonal) will be zero
zero_ind = (index_i == index_j-1) | (index_i == index_j+1);
index_i(zero_ind) = NaN;
index_j(zero_ind) = NaN;

% indexing masks for upper and lower triangular part
data_L_mask = index_i >= index_j;
data_U_mask = index_i < index_j;

data_L_i = index_i(data_L_mask);
data_L_j = index_j(data_L_mask);
data_U_i = index_i(data_U_mask);
data_U_j = index_j(data_U_mask);

% extract the linear indeces from the data matrix and regularization matrix
WW_L_ind = find(WW_L);
data_L_ind = sub2ind(size(WW_L), data_L_i, data_L_j);

% compute the linear indeces of the non-zeros in the matrix
[L_ind, ~, data_WW_in_L_index] = unique([data_L_ind; WW_L_ind]);

% compute the corresponding indices for the values in the data and reg
% matrix
data_in_L_index = uint32(data_WW_in_L_index(1:length(data_L_ind)));
WW_in_L_index = data_WW_in_L_index(length(data_L_ind)+1:end);

% create the arrays of values in the regularization matrix
nnz_L = length(L_ind);
WW_L_vec = zeros(nnz_L, 1, 'single');
WW_L_vec(WW_in_L_index) = full(WW_L(WW_L_ind));

% precompute the data part of the regularization matrix
WW_L_vec_data = WW_L_vec(data_in_L_index);

% initialize the content vectors for the matrices
L_vec = WW_L_vec;

% preallocate the matrices
mat_size = feature_dim * support_sz;
[L_i, L_j] = ind2sub(size(WW_L), L_ind);
AL = sparse(L_i, L_j, ones(nnz_L,1), mat_size, mat_size);
AU_data = sparse(data_U_i, data_U_j, ones(length(data_U_i),1), mat_size, mat_size);

if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end

% initialize the projection matrix
rect_position = zeros(num_frames, 4);

time = 0;

% allocate
xxlf_sep = zeros(2*feature_dim, length(dft_sym_ind) + 2 * length(dft_pos_ind), feature_dim, 'single');
multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');

%pos_diff = zeros(num_frames,2);
key_1 = 0;
key_2 = 0;
key_3 = 0;
threshold =12.38;
test_t = zeros(num_frames,6);
for frame = 1:num_frames,
    %load image
    im = imread(s_frames{frame});
    if size(im,3) > 1 && colorImage == false
        im = im(:,:,1);
    end
	learning_rate_(1) = learning_rate;
	learning_rate_(2) = learning_rate;
	learning_rate_(3) = learning_rate;
    tic();
    
    %do not estimate translation and scaling on the first frame, since we 
    %just want to initialize the tracker there
    yf_vec1 = yf_vec_1;
    yf_vec2 = yf_vec_2;
    yf_vec3 = yf_vec_3;
    if frame > 1
        old_pos = inf(size(pos));
        iter = 1;
        
        %translation search
        while iter <= refinement_iterations && any(old_pos ~= pos)
            % Get multi-resolution image
            for scale_ind = 1:nScales
                multires_pixel_template(:,:,:,scale_ind) = ...
                    get_pixels(im, pos, round(sz*currentScaleFactor*scaleFactors(scale_ind)), sz);
            end
            
            xt = bsxfun(@times,get_features(multires_pixel_template,features,global_feat_params),cos_window);
            
            xtf = fft2(xt);

            responsef1 = permute(sum(bsxfun(@times, hf1, xtf), 3), [1 2 4 3]);
            responsef2 = permute(sum(bsxfun(@times, hf2, xtf), 3), [1 2 4 3]);
            responsef3 = permute(sum(bsxfun(@times, hf3, xtf), 3), [1 2 4 3]);
            % if we undersampled features, we want to interpolate the
            % response so it has the same size as the image patch
            if interpolate_response == 2
                % use dynamic interp size
                interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
            end
            responsef_padded1 = resizeDFT2(responsef1, interp_sz);
            responsef_padded2 = resizeDFT2(responsef2, interp_sz);
            responsef_padded3 = resizeDFT2(responsef3, interp_sz);

            % response
            response1 = ifft2(responsef_padded1, 'symmetric');
            response2 = ifft2(responsef_padded2, 'symmetric');
            response3 = ifft2(responsef_padded3, 'symmetric');

            % find maximum
            if interpolate_response == 3
                error('Invalid parameter value for interpolate_response');
            elseif interpolate_response == 4
                [disp_row1, disp_col1, sind1,max_scale_response1] = resp_newton(response1, responsef_padded1, newton_iterations, ky, kx, use_sz);
				[disp_row2, disp_col2, sind2,max_scale_response2] = resp_newton(response2, responsef_padded2, newton_iterations, ky, kx, use_sz);
                [disp_row3, disp_col3, sind3,max_scale_response3] = resp_newton(response3, responsef_padded3, newton_iterations, ky, kx, use_sz);
            %get tran
                disp_row1 = disp_row1+1;disp_col1 = disp_col1-3;
                disp_row2 = disp_row2-1;disp_col2 = disp_col2+3;
            else
                [row1, col1, sind1] = ind2sub(size(response1), find(response1 == max(response1(:)), 1));
				[row2, col2, sind2] = ind2sub(size(response2), find(response2 == max(response2(:)), 1));
                [row3, col3, sind3] = ind2sub(size(response3), find(response3 == max(response3(:)), 1));
                disp_row1 = mod(row1 - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
				disp_row2 = mod(row2 - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
                disp_row3 = mod(row3 - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
                disp_col1 = mod(col1 - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
				disp_col2 = mod(col2 - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
                disp_col3 = mod(col3 - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
            end

            % calculate translation
            switch interpolate_response
                case 0
                    translation_vec1 = round([disp_row1, disp_col1] * featureRatio * currentScaleFactor * scaleFactors(sind1));
					            subS = subSeqs{1};
                                subSeqs=[];
                                subSeqs{1} = subS;

                                subA = subAnno{1};
                                subAnno=[];
                                subAnno{1} = subA;
                case 1
                    translation_vec1 = round([disp_row1, disp_col1] * currentScaleFactor * scaleFactors(sind1));
					translation_vec2 = round([disp_row2, disp_col2] * currentScaleFactor * scaleFactors(sind2));
					translation_vec3 = round([disp_row3, disp_col3] * currentScaleFactor * scaleFactors(sind2));
                case 2
                    translation_vec1 = round([disp_row1, disp_col1] * scaleFactors(sind1));
					translation_vec2 = round([disp_row2, disp_col2] * scaleFactors(sind2));
					translation_vec3 = round([disp_row3, disp_col3] * scaleFactors(sind3));
                case 3
                    translation_vec1 = round([disp_row1, disp_col1] * featureRatio * currentScaleFactor * scaleFactors(sind1));
					translation_vec2 = round([disp_row2, disp_col2] * featureRatio * currentScaleFactor * scaleFactors(sind2));
					translation_vec3 = round([disp_row3, disp_col3] * featureRatio * currentScaleFactor * scaleFactors(sind3));
                case 4
                    translation_vec1 = round([disp_row1, disp_col1] * featureRatio * currentScaleFactor * scaleFactors(sind1));
					translation_vec2 = round([disp_row2, disp_col2] * featureRatio * currentScaleFactor * scaleFactors(sind2));
					translation_vec3 = round([disp_row3, disp_col3] * featureRatio * currentScaleFactor * scaleFactors(sind3));
            end

		%##evaluate
			%##evaluate PSR
				scale_response1 = response1(:,:,sind1);
				scale_response2 = response2(:,:,sind2);
				scale_response3 = response3(:,:,sind3);
				
				average_scale_response1 = mean(scale_response1(:));
				average_scale_response2 = mean(scale_response2(:));
				average_scale_response3 = mean(scale_response3(:));
				
				variance_response1 = std2(scale_response1);
				variance_response2 = std2(scale_response2);
				variance_response3 = std2(scale_response3);
				
				PSR1 = (max_scale_response1-average_scale_response1)/variance_response1;
				PSR2 = (max_scale_response2-average_scale_response2)/variance_response2;
				PSR3 = (max_scale_response3-average_scale_response3)/variance_response3;
			%##end evaluate PSR

				%##conserve response to old  and  row/col to old
					old_scale_response1 = scale_response1;
					old_scale_response2 = scale_response2;
					old_scale_response3 = scale_response3;
				%## end conserve
				weight1 = PSR1;
				weight2 = PSR2;
				weight3 = PSR3;
			%##end combine to weight
			%##mark of no occlude
			mark1=1;mark2=1;mark3=1;
 			%max_weight = max([weight1 weight2 weight3]);
				if weight1<threshold
					%mark1 = 0;
                    learning_rate_(1)=learning_rate_(1)/2;
                    max_scale_response1 = max_scale_response1/2;
				end
				if weight2<threshold
					%mark2 = 0;
					learning_rate_(2)=learning_rate_(2)/2;
                    max_scale_response2 = max_scale_response2/2;
				end
				if weight3<threshold
					%mark3 = 0;
					learning_rate_(3)=learning_rate_(3)/2;
                    max_scale_response3 = max_scale_response3/2;
				end
			%##end mark of occlude
		%##end evaluate
		
		%##scale evaluate 
		

			sind_0 = old_sind_0;
			if (mark1+mark2+mark3)~=0
				sind_0 = round((mark1*sind1 + mark2*sind2 + mark3*sind3)/(mark1+mark2+mark3));
			end
			old_sind_0 = sind_0;
		%##end scale evaluate
		
            % set the scale
            currentScaleFactor = currentScaleFactor * scaleFactors(sind_0);
			
            % adjust to make sure we are not to large or to small
            if currentScaleFactor < min_scale_factor
                currentScaleFactor = min_scale_factor;
            elseif currentScaleFactor > max_scale_factor
                currentScaleFactor = max_scale_factor;
            end
	
%             end
			target_sz = floor(base_target_sz * currentScaleFactor);%	pos-sz+translation<0
			test = pos() - floor(target_sz/2);
			test_object1 = test + translation_vec1;
			test_object2 = test + translation_vec2;
			test_object3 = test + translation_vec3;
			for i = 1:2
				if test_object1(i)<0
					translation_vec1(i) = -test(i);
					learning_rate_(1) = 0;
                   % disp(frame);
				end
				if test_object2(i)<0
					translation_vec2(i) = -test(i);
					learning_rate_(2) = 0;
                    %disp(frame);
				end
				if test_object3(i)<0
					translation_vec3(i) = -test(i);
					learning_rate_(3) = 0;
                    %disp(frame);
				end
			end
		
		
			trans = 0;
			if (mark1+mark2+mark3)~=0
				trans = (max_scale_response1*mark1*translation_vec1 + ...
							max_scale_response2*mark2*translation_vec2 + ...
								max_scale_response3*mark3*translation_vec3)/...
						(max_scale_response1*mark1 + max_scale_response2*mark2 + max_scale_response3*mark3);
			end
			pos = pos + round(trans);
            iter = iter + 1;
        end

        % debug visualization of responses
        if debug
            figure(101);
            subplot_cols = ceil(sqrt(nScales));
            subplot_rows = ceil(nScales/subplot_cols);
            for scale_ind = 1:nScales
                subplot(subplot_rows,subplot_cols,scale_ind);
                imagesc(fftshift(response(:,:,scale_ind)));colorbar; axis image;
                title(sprintf('Scale %i,  max(response) = %f', scale_ind, max(max(response(:,:,scale_ind)))));
            end
        end
    end
    
    % extract training sample image region
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    
    % extract features and do windowing
    xl = bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window);
    
    % take the DFT and vectorize each feature dimension
    xlf = fft2(xl);
    xlf_reshaped = reshape(xlf, [support_sz, feature_dim]);
    
	
    % new rhs sample
    xyf_corr1 = bsxfun(@times, yf_vec1, conj(xlf_reshaped));
	xyf_corr2 = bsxfun(@times, yf_vec2, conj(xlf_reshaped));
	xyf_corr3 = bsxfun(@times, yf_vec3, conj(xlf_reshaped));
    xy_dfs1 = real(dfs_matrix * double(xyf_corr1));
	xy_dfs2 = real(dfs_matrix * double(xyf_corr2));
	xy_dfs3 = real(dfs_matrix * double(xyf_corr3));
	
    new_hf_rhs1 = xy_dfs1(:);
    new_hf_rhs2 = xy_dfs2(:);
	new_hf_rhs3 = xy_dfs3(:);
	
    xlf_reshaped_sym = xlf_reshaped(dft_sym_ind, :);    % extract the symmetric part of the spectrum x_0
    xlf_reshaped_pos = xlf_reshaped(dft_pos_ind, :);    % extract the positive part of the spectrum x_+
    
    % compute autocorrelation
    xxlf_sym = bsxfun(@times, conj(permute(xlf_reshaped_sym, [2 1])), permute(xlf_reshaped_sym, [3 1 2]));
    xxlf_pos = bsxfun(@times, conj(permute(xlf_reshaped_pos, [2 1])), permute(xlf_reshaped_pos, [3 1 2]));
    xxlf_pos_real = real(xxlf_pos);
    
    % partition the real and imaginary parts
    xxlf_sep(1:2:end, 1:num_sym_coef, :) = real(xxlf_sym);
    xxlf_sep(1:2:end, num_sym_coef+1:2:end, :) = xxlf_pos_real;
    xxlf_sep(2:2:end, num_sym_coef+1:2:end, :) = imag(xxlf_pos);
    xxlf_sep(1:2:end, num_sym_coef+2:2:end, :) = -imag(xxlf_pos);
    xxlf_sep(2:2:end, num_sym_coef+2:2:end, :) = xxlf_pos_real;
    
    if frame == 1
        hf_rhs1 = new_hf_rhs1;
		hf_rhs2 = new_hf_rhs2;
		hf_rhs3 = new_hf_rhs3;
        hf_autocorr = xxlf_sep(:);
        
        % compute the initial filter in the first frame
        hf_init_autocorr = double(sum(xlf_reshaped .* conj(xlf_reshaped), 2));
        switch params.init_strategy
            case 'const_reg'       % exact solution for constant regularization
                hf_init = bsxfun(@rdivide, xyf_corr, hf_init_autocorr + params.reg_window_min^2);
                hf_init = real(dfs_matrix * hf_init);
                hf_vec = hf_init(:);
            case 'indep'           % independent filters for each feature
                A_init = real(dfs_matrix * spdiags(hf_init_autocorr, 0, support_sz, support_sz) * dfs_matrix') + feature_dim * WW_block;
                b_init1 = reshape(hf_rhs1, support_sz, feature_dim);
				b_init2 = reshape(hf_rhs2, support_sz, feature_dim);
				b_init3 = reshape(hf_rhs3, support_sz, feature_dim);
                hf_init1 = A_init \ b_init1;
				hf_init2 = A_init \ b_init2;
				hf_init3 = A_init \ b_init3;
                hf_vec1 = hf_init1(:);
				hf_vec2 = hf_init2(:);
				hf_vec3 = hf_init3(:);
        end
		%##update old information 
			old_scale_response1 = y1;
			old_scale_response2 = y2;
			old_scale_response3 = y3;
			old_sind_0 = 4;
		%## end
        
        if debug == 1
            hf_vec_old = hf_vec;
            hf_difference = zeros(num_frames * num_GS_iter, 1);
        end
    else
        hf_rhs1 = (1 - learning_rate_(1)) * hf_rhs1 + learning_rate_(1) * new_hf_rhs1;				
		hf_rhs2 = (1 - learning_rate_(2)) * hf_rhs2 + learning_rate_(2) * new_hf_rhs2;				
		hf_rhs3 = (1 - learning_rate_(3)) * hf_rhs3 + learning_rate_(3) * new_hf_rhs3;
        hf_autocorr = (1 - learning_rate) * hf_autocorr + learning_rate * xxlf_sep(:);
		
    end
    
    % add the autocorrelation to the matrix vectors with the regularization
    L_vec(data_in_L_index) = hf_autocorr(data_L_mask) + WW_L_vec_data;
    
    % update the matrices with the new non-zeros
    AL = setnonzeros(AL, double(L_vec));
    AU_data = setnonzeros(AU_data, double(hf_autocorr(data_U_mask)));
    
    % do Gausss-Seidel
    for iter = 1:num_GS_iter
        hf_vec1 = AL \ (hf_rhs1 - AU_data * hf_vec1 - WW_U * hf_vec1);
        hf_vec2 = AL \ (hf_rhs2 - AU_data * hf_vec2 - WW_U * hf_vec2);
		hf_vec3 = AL \ (hf_rhs3 - AU_data * hf_vec3 - WW_U * hf_vec3);
        if debug
            hf_difference((frame - 1)*num_GS_iter + iter) = norm(hf_vec - hf_vec_old);
            hf_vec_old1 = hf_vec1;
			hf_vec_old2 = hf_vec2;
			hf_vec_old3 = hf_vec3;
        end
    end
    
    % reconstruct the filter
    hf1 = reshape(single(dfs_matrix' * reshape(hf_vec1, [support_sz, feature_dim])), [use_sz, feature_dim]);
    hf2 = reshape(single(dfs_matrix' * reshape(hf_vec2, [support_sz, feature_dim])), [use_sz, feature_dim]);
	hf3 = reshape(single(dfs_matrix' * reshape(hf_vec3, [support_sz, feature_dim])), [use_sz, feature_dim]);

   % debug visualization
    if debug
        % plot filters
        figure(20)
        subplot_cols = ceil(sqrt(feature_dim));
        subplot_rows = ceil(feature_dim/subplot_cols);
        for disp_layer = 1:feature_dim
            subplot(subplot_rows,subplot_cols,disp_layer);
            imagesc(ifft2(conj(hf(:,:,disp_layer)), 'symmetric')); 
            colorbar;
            axis image;
        end
        
        % plot convergence
        figure(99);plot(hf_difference);axis([1, num_GS_iter * frame, 0, max(hf_difference)]);
        title('Filter optimization convergence');
    end
    
    target_sz = floor(base_target_sz * currentScaleFactor);
    
    %save position and calculate FPS
    rect_position(frame,:) = [pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
	time = time + toc();
    
    %visualization
    if visualization == 1
        rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        im_to_show = double(im)/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == 1
            fig_handle = figure('Name', 'Tracking');
            imagesc(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
            axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
        else
            resp_sz = round(sz*currentScaleFactor*scaleFactors(scale_ind));
            xs = floor(old_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
            ys = floor(old_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
            sc_ind = floor((nScales - 1)/2) + 1;
            
            figure(fig_handle);
            imagesc(im_to_show);
            hold on;
            resp_handle1 = imagesc(xs, ys, fftshift(response1(:,:,sc_ind))); colormap hsv;
            resp_handle2= imagesc(xs, ys, fftshift(response2(:,:,sc_ind))); colormap hsv;
            resp_handle3 = imagesc(xs, ys, fftshift(response3(:,:,sc_ind))); colormap hsv;
            alpha(resp_handle1, 0.25);
            alpha(resp_handle2, 0.25);
            alpha(resp_handle3, 0.25);
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
        end
        
        drawnow
        % pause
    end
end

fps = numel(s_frames) / time;

results.type = 'rect';
results.res = rect_position;
results.fps = fps;
results.test = test_t;
%results.pos_diff = pos_diff;