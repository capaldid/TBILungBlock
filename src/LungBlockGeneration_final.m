%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:
%   This function performs the automatic generation of lung block shells
%   to be filled with tungsten ball bearings from radiation therapy
%   treatment plan DICOM files
%
% Inputs:
%         RT Plan dcm file
% 
% Outputs:
%         STL files for all four lung blocks shells (tops and bottoms)
% 
% Author: Dante PI Capaldi
% Website: https://github.com/capaldid
% Date: July 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


close all; clear all; clc;

% parameters (information is in mm)
shift_org = 200; % shift is used to create the stl file
height_block = 15+6; % height of lung block (mm)
screw_dia = 5; % screw  diameter
wall_thickness = 2; % thickness of wall
slit_separation = 32; % distance between upper and lower slit to hold blocks
iso_distance = 30; % max distance between iso marking and center of screw positions
min_size_block = 50; % minimum size of lung blocks (was previously 70 mm, had to change for patient 11 on July 13 2023)
min_distance_edge = 15; % minimum distance between screw holes and outer shell


% select RT Plan File
[file,path] = uigetfile('*.dcm');
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

% running the stl code
if isequal(file,0)
    % do not run code and cancel software
    disp('Code Cancel');
else
    % read dcm header
    H = dicominfo([path file]);
    
    % extracting first lung block trace from dcm header
    tmp = H.BeamSequence.Item_1.BlockSequence.Item_1.BlockData;
    x = round(tmp(1:2:(size(tmp,1)-1)))+shift_org;
    y = round(tmp(2:2:(size(tmp,1))))+shift_org;
    
    % extracting second lung block trace from dcm header
    tmp2 = H.BeamSequence.Item_1.BlockSequence.Item_2.BlockData;
    x2 = round(tmp2(1:2:(size(tmp2,1)-1)))+shift_org;
    y2 = round(tmp2(2:2:(size(tmp2,1))))+shift_org;
    
    % plot first set of lung blocks
    figure,
    subplot(1,2,1)
    plot(x,y,x2,y2);
    title(H.BeamSequence.Item_1.BeamName)
    
    % extracting third lung block trace from dcm header
    tmp3 = H.BeamSequence.Item_2.BlockSequence.Item_1.BlockData;
    x3 = round(tmp3(1:2:(size(tmp3,1)-1)))+shift_org;
    y3 = round(tmp3(2:2:(size(tmp3,1))))+shift_org;
    
    % extracting fourth lung block trace from dcm header
    tmp4 = H.BeamSequence.Item_2.BlockSequence.Item_2.BlockData;
    x4 = round(tmp4(1:2:(size(tmp4,1)-1)))+shift_org;
    y4 = round(tmp4(2:2:(size(tmp4,1))))+shift_org;
    
    % plot second set of lung blocks
    subplot(1,2,2)
    plot(x3,y3,x4,y4);
    title(H.BeamSequence.Item_2.BeamName)

    % producing contours from lung block traces
    % note: to calc volume of shell = sum(BW(:))*19/1000
    % 19 is the height of shell and 1000 is to convert to cc
    tmpImg = zeros(shift_org*2,shift_org*2);
    BW1 = roipoly(tmpImg,x,y);
    BW2 = roipoly(tmpImg,x2,y2);
    BW3 = roipoly(tmpImg,x3,y3);
    BW4 = roipoly(tmpImg,x4,y4);

    % get stats of the centroids
    BW1_stats = regionprops(BW1);
    BW1_centroid = round(BW1_stats.Centroid);
    BW2_stats = regionprops(BW2);
    BW2_centroid = round(BW2_stats.Centroid);
    BW3_stats = regionprops(BW3);
    BW3_centroid = round(BW3_stats.Centroid);
    BW4_stats = regionprops(BW4);
    BW4_centroid = round(BW4_stats.Centroid);
    
    % dilate the borders on the thickness of the walls
    se = strel('disk',wall_thickness);
    BW1_dil = double(imdilate(BW1,se));
    BW2_dil = double(imdilate(BW2,se));
    BW3_dil = double(imdilate(BW3,se));
    BW4_dil = double(imdilate(BW4,se));

    % determining the distance between the center of the screw holes and
    % the edges (sup/inf direction) to make sure the screw holes are not
    % outside of the block (for the first set of blocks)

    % for the supine blocks

    % ... for the first lung block, sup direction
    tmp_BW1 = BW1_dil-double(BW1);
    flag = 0; n1 = 0;
    while flag == 0
        if tmp_BW1(BW1_centroid(2)+n1,BW1_centroid(1)) == 0
            n1 = n1 + 1;
        else
            flag = 1;
        end
    end
    % ... for the first lung block, inf direction
    n2 = 0; flag = 0;
    while flag == 0
        if tmp_BW1(BW1_centroid(2)-n2,BW1_centroid(1)) == 0
            n2 = n2 + 1;
        else
            flag = 1;
        end
    end
    % determine the distances above/below centroid for first lung block
    BW1_bottom_val = BW1_centroid(2)-n2;
    BW1_top_val = BW1_centroid(2)+n1;

    % ... for the second lung block, sup direction
    tmp_BW2 = BW2_dil-double(BW2);
    flag = 0; n1 = 0;
    while flag == 0
        if tmp_BW2(BW2_centroid(2)+n1,BW2_centroid(1)) == 0
            n1 = n1 + 1;
        else
            flag = 1;
        end
    end
    % ... for the second lung block, inf direction
    flag = 0; n2 = 0;
    while flag == 0
        if tmp_BW2(BW2_centroid(2)-n2,BW2_centroid(1)) == 0
            n2 = n2 + 1;
        else
            flag = 1;
        end
    end

    % determine the distances above/below centroid for second lung block
    BW2_bottom_val = BW2_centroid(2)-n2;
    BW2_top_val = BW2_centroid(2)+n1;

    % finding the minimum distances for the lung block calc (supine)
    [val1,idx1]= min([BW1_top_val BW2_top_val]);
    [val2,idx2]= max([BW1_bottom_val BW2_bottom_val]);

    % for the prone blocks

    % ... for the first lung block, sup direction
    tmp_BW3 = BW3_dil-double(BW3);
    flag = 0; n3 = 0;
    while flag == 0
        if tmp_BW3(BW3_centroid(2)+n3,BW3_centroid(1)) == 0
            n3 = n3 + 1;
        else
            flag = 1;
        end
    end
    % ... for the first lung block, inf direction
    n4 = 0; flag = 0;
    while flag == 0
        if tmp_BW3(BW3_centroid(2)-n4,BW3_centroid(1)) == 0
            n4 = n4 + 1;
        else
            flag = 1;
        end
    end
    % determine the distances above/below centroid for first lung block
    BW3_bottom_val = BW3_centroid(2)-n4;
    BW3_top_val = BW3_centroid(2)+n3;

    % ... for the second lung block, sup direction
    tmp_BW4 = BW4_dil-double(BW4);
    flag = 0; n3 = 0;
    while flag == 0
        if tmp_BW4(BW4_centroid(2)+n3,BW4_centroid(1)) == 0
            n3 = n3 + 1;
        else
            flag = 1;
        end
    end
    % ... for the second lung block, inf direction
    flag = 0; n4 = 0;
    while flag == 0
        if tmp_BW4(BW4_centroid(2)-n4,BW4_centroid(1)) == 0
            n4 = n4 + 1;
        else
            flag = 1;
        end
    end

    % determine the distances above/below centroid for second lung block
    BW4_bottom_val = BW4_centroid(2)-n4;
    BW4_top_val = BW4_centroid(2)+n3;

    % finding the minimum distances for the lung block calc (prone)
    [val3,idx3]= min([BW3_top_val BW4_top_val]);
    [val4,idx4]= max([BW3_bottom_val BW4_bottom_val]);

    % if the lung blocks are too small, we cannot print them
    if (val1 - val2) < min_size_block || (val3 - val4) < min_size_block
        disp("Cannot create lung blocks, too small to create")
    else
        % Circle code to create the screw cutout
        r1 = round(screw_dia/2);
        x = 1:size(tmpImg,2);
        y = 1:size(tmpImg,1);
        [xx,yy] = meshgrid(x,y);
    
        % shift value is a parameter used to shift where the screw holes
        % are relative to the lung blocks, based on the geometry of the
        % blocks
        shift_val_supine = 0;

        % based on the geometry, the shift value will change to satisfy all
        % criteria
        for i = 0:100
            
            % find the average centroid location between the left and right
            % lung blocks
            avg_BW1BW2_centroid = round((BW1_centroid(2)+BW2_centroid(2))/2)+ shift_val_supine;
            xc_1 = BW1_centroid(1);
            xc_2 = BW2_centroid(1);
        
            % shift the centroid locations up and down to where the screw
            % locations should be (left lung block)
            yc_1a = avg_BW1BW2_centroid-round(slit_separation/2);
            yc_1b = avg_BW1BW2_centroid+round(slit_separation/2);

            % shift the centroid locations up and down to where the screw
            % locations should be (right lung block)
            yc_2a = avg_BW1BW2_centroid-round(slit_separation/2);
            yc_2b = avg_BW1BW2_centroid+round(slit_separation/2);

            
            % determine if there is enough room on the transparent bridge
            % between the isocenter on the blocks and the screw locations
            % on the block
            if abs(avg_BW1BW2_centroid - shift_org) < iso_distance
                % determine of the screw holes are appropriately placed on
                % the lung blocks based on the criteria of the
                % min_distance_edge (both left and right lung blocks, sup
                % and inf direction).  If it doesn't, shift the screw
                % locations so that they are satisfied
                if (yc_1a - BW1_bottom_val) < 0 || abs(yc_1a - BW1_bottom_val) < min_distance_edge
                    shift_val_supine = shift_val_supine + 1;
                elseif (yc_1b - BW1_top_val) > 0 || abs(yc_1b - BW1_top_val) < min_distance_edge
                    shift_val_supine = shift_val_supine - 1;
                elseif (yc_2a - BW2_bottom_val) < 0 || abs(yc_2a - BW2_bottom_val) < min_distance_edge
                    shift_val_supine = shift_val_supine + 1;
                elseif (yc_2b - BW2_top_val) > 0 || abs(yc_2b - BW2_top_val) < min_distance_edge
                    shift_val_supine = shift_val_supine - 1;
                else
                    % if all criteria is satisfied, draw screw locations on
                    % left lung block
                    mask1_1a = hypot(xx - xc_1, yy - yc_1a) <= r1;
                    mask1_1b = hypot(xx - xc_1, yy - yc_1b) <= r1;
                    BW1_total = (BW1_dil-double(BW1));
                    % if all criteria is satisfied, draw screw locations on
                    % right lung block
                    mask1_2a = hypot(xx - xc_2, yy - yc_2a) <= r1;
                    mask1_2b = hypot(xx - xc_2, yy - yc_2b) <= r1;
                    BW2_total = (BW2_dil-double(BW2));
                    break
                end
        
            else
                % if it doesn't satisfy isocenter criteria, perform shift
                if (avg_BW1BW2_centroid - shift_org) < 0
                    shift_val_supine = shift_val_supine + 1;
                else
                    shift_val_supine = shift_val_supine - 1;
                end
            end

        end

        % after 100 iterations, if not satisfied, find the optimal solution
        if i == 100
            mask1_1a = hypot(xx - xc_1, yy - yc_1a) <= r1;
            mask1_1b = hypot(xx - xc_1, yy - yc_1b) <= r1;
            BW1_total = (BW1_dil-double(BW1));

            mask1_2a = hypot(xx - xc_2, yy - yc_2a) <= r1;
            mask1_2b = hypot(xx - xc_2, yy - yc_2b) <= r1;
            BW2_total = (BW2_dil-double(BW2));
        end

        % display final masks
        figure,
        subplot(3,2,1.5)
        imshow(rot90(rot90(BW2_total+BW1_total)),[])
        title(H.BeamSequence.Item_1.BeamName)

        % perform dilation for block tops
        se = strel('disk',1);
        BW1_dil1 = double(imdilate(BW1,se));
        BW2_dil1 = double(imdilate(BW2,se));

        % generate block tops based on dilation (and erode for fitting
        % purposes)
        BW_sup_vol_top = zeros(shift_org*2,shift_org*2,5);
        BW_sup_vol_top(:,:,2) = imerode(BW1 + BW2,se) - (mask1_1a + mask1_1b) - (mask1_2a + mask1_2b);
        BW_sup_vol_top(:,:,3) = imerode(BW1_dil1 + BW2_dil1,se)  - (mask1_1a + mask1_1b) - (mask1_2a + mask1_2b);
        BW_sup_vol_top(:,:,4) = imerode(BW1 + BW2,se) - (mask1_1a + mask1_1b) - (mask1_2a + mask1_2b);
    
        % producing supine lung block volume
        BW_sup_vol = repmat((BW1_total+BW2_total),1,1,height_block+5);
        BW_sup_vol(:,:,1) = 0; % remove bottom slice
        % generate the bottom for the lung block with screw holes
        for i = 1:wall_thickness
            BW_sup_vol(:,:,i+1) = BW2_dil + BW1_dil - mask1_1a - mask1_2a - mask1_1b - mask1_2b;
        end

        % create notch on the top of the bottom block so that top fits in
        BW_sup_vol(:,:,height_block+3) = BW2_dil - double(BW2) + BW1_dil-double(BW1) - (BW1_dil1 + BW2_dil1);
        BW_sup_vol(:,:,height_block+4) = BW2_dil - double(BW2) + BW1_dil-double(BW1) - (BW1 + BW2);
        BW_sup_vol(:,:,height_block+5) = 0;

        % tmp code to do only half block
        BW_sup_vol_left = BW_sup_vol;
        BW_sup_vol_top_left = BW_sup_vol_top;
        BW_sup_vol_left(:,201:end,:) = 0;
        BW_sup_vol_top_left(:,201:end,:) = 0;

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_sup_vol_left,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_sup_vol_left,1),size(BW_sup_vol_left,2),size(BW_sup_vol_left,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        % % code to display 3D rendering appropriately (not used here)
        % viewerIntrinsic = viewer3d;
        % obj = images.ui.graphics3d.Surface(viewerIntrinsic, ...
        % Data=triInstrinsic, ...
        % Color=[0.88 0.84 0.71], ...
        % Alpha=0.9);

        subplot(3,2,3)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.1])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_1.BeamName 'BottomLeft.stl'])

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_sup_vol_top_left,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_sup_vol_top_left,1),size(BW_sup_vol_top_left,2),size(BW_sup_vol_top_left,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,4)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.01])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_1.BeamName 'TopLeft.stl'])

        % tmp code to do only half block
        BW_sup_vol_right = BW_sup_vol;
        BW_sup_vol_top_right = BW_sup_vol_top;
        BW_sup_vol_right(:,1:200,:) = 0;
        BW_sup_vol_top_right(:,1:200,:) = 0;

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_sup_vol_right,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_sup_vol_right,1),size(BW_sup_vol_right,2),size(BW_sup_vol_right,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,5)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.1])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_1.BeamName 'BottomRight.stl'])

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_sup_vol_top_right,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_sup_vol_top_right,1),size(BW_sup_vol_top_right,2),size(BW_sup_vol_top_right,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,6)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.01])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_1.BeamName 'TopRight.stl'])


        % shift value is a parameter used to shift where the screw holes
        % are relative to the lung blocks, based on the geometry of the
        % blocks
        shift_val_prone = 0;

        % based on the geometry, the shift value will change to satisfy all
        % criteria
        for i = 0:100
            
            % find the average centroid location between the left and right
            % lung blocks
            avg_BW3BW4_centroid = round((BW3_centroid(2)+BW4_centroid(2))/2)+ shift_val_prone;
            xc_3 = BW3_centroid(1);
            xc_4 = BW4_centroid(1);
        
            % shift the centroid locations up and down to where the screw
            % locations should be (left lung block)
            yc_3a = avg_BW3BW4_centroid-round(slit_separation/2);
            yc_3b = avg_BW3BW4_centroid+round(slit_separation/2);

            % shift the centroid locations up and down to where the screw
            % locations should be (right lung block)
            yc_4a = avg_BW3BW4_centroid-round(slit_separation/2);
            yc_4b = avg_BW3BW4_centroid+round(slit_separation/2);

            
            % determine if there is enough room on the transparent bridge
            % between the isocenter on the blocks and the screw locations
            % on the block
            if abs(avg_BW3BW4_centroid - shift_org) < iso_distance
                % determine of the screw holes are appropriately placed on
                % the lung blocks based on the criteria of the
                % min_distance_edge (both left and right lung blocks, sup
                % and inf direction).  If it doesn't, shift the screw
                % locations so that they are satisfied
                if (yc_3a - BW3_bottom_val) < 0 || abs(yc_3a - BW3_bottom_val) < min_distance_edge
                    shift_val_prone = shift_val_prone + 1;
                elseif (yc_3b - BW3_top_val) > 0 || abs(yc_3b - BW3_top_val) < min_distance_edge
                    shift_val_prone = shift_val_prone - 1;
                elseif (yc_4a - BW4_bottom_val) < 0 || abs(yc_4a - BW4_bottom_val) < min_distance_edge
                    shift_val_prone = shift_val_prone + 1;
                elseif (yc_4b - BW4_top_val) > 0 || abs(yc_4b - BW4_top_val) < min_distance_edge
                    shift_val_prone = shift_val_prone - 1;
                else
                    % if all criteria is satisfied, draw screw locations on
                    % left lung block
                    mask3_1a = hypot(xx - xc_3, yy - yc_3a) <= r1;
                    mask3_1b = hypot(xx - xc_3, yy - yc_3b) <= r1;
                    BW3_total = (BW3_dil-double(BW3));
                    % if all criteria is satisfied, draw screw locations on
                    % right lung block
                    mask3_2a = hypot(xx - xc_4, yy - yc_4a) <= r1;
                    mask3_2b = hypot(xx - xc_4, yy - yc_4b) <= r1;
                    BW4_total = (BW4_dil-double(BW4));
                    break
                end
        
            else
                % if it doesn't satisfy isocenter criteria, perform shift
                if (avg_BW3BW4_centroid - shift_org) < 0
                    shift_val_prone = shift_val_prone + 1;
                else
                    shift_val_prone = shift_val_prone - 1;
                end
            end

        end

        % after 100 iterations, if not satisfied, find the optimal solution
        if i == 100
            mask3_1a = hypot(xx - xc_3, yy - yc_3a) <= r1;
            mask3_1b = hypot(xx - xc_3, yy - yc_3b) <= r1;
            BW3_total = (BW3_dil-double(BW3));

            mask3_2a = hypot(xx - xc_4, yy - yc_4a) <= r1;
            mask3_2b = hypot(xx - xc_4, yy - yc_4b) <= r1;
            BW4_total = (BW4_dil-double(BW4));
        end

        % display final masks
        figure,
        subplot(3,2,1.5)
        imshow(rot90(rot90(BW3_total+BW4_total)),[])
        title(H.BeamSequence.Item_2.BeamName)

        % perform dilation for block tops
        se = strel('disk',1);
        BW3_dil1 = double(imdilate(BW3,se));
        BW4_dil1 = double(imdilate(BW4,se));

        % generate block tops based on dilation (and erode for fitting
        % purposes)
        BW_pro_vol_top = zeros(shift_org*2,shift_org*2,5);
        BW_pro_vol_top(:,:,2) = imerode(BW3 + BW4,se) - (mask3_1a + mask3_1b) - (mask3_2a + mask3_2b);
        BW_pro_vol_top(:,:,3) = imerode(BW3_dil1 + BW4_dil1,se) - (mask3_1a + mask3_1b) - (mask3_2a + mask3_2b);
        BW_pro_vol_top(:,:,4) = imerode(BW3 + BW4,se) - (mask3_1a + mask3_1b) - (mask3_2a + mask3_2b);
    
        % producing supine lung block volume
        BW_pro_vol = repmat((BW3_total+BW4_total),1,1,height_block+5);
        BW_pro_vol(:,:,1) = 0; % remove bottom slice
        % generate the bottom for the lung block with screw holes
        for i = 1:wall_thickness
            BW_pro_vol(:,:,i+1) = BW4_dil + BW3_dil - mask3_1a - mask3_2a - mask3_1b - mask3_2b;
        end

        % create notch on the top of the bottom block so that top fits in
        BW_pro_vol(:,:,height_block+3) = BW4_dil - double(BW4) + BW3_dil-double(BW3) - (BW3_dil1 + BW4_dil1);
        BW_pro_vol(:,:,height_block+4) = BW4_dil - double(BW4) + BW3_dil-double(BW3) - (BW3 + BW4);
        BW_pro_vol(:,:,height_block+5) = 0;

        % tmp code to do only half block
        BW_pro_vol_left = BW_pro_vol;
        BW_pro_vol_top_left = BW_pro_vol_top;
        BW_pro_vol_left(:,201:end,:) = 0;
        BW_pro_vol_top_left(:,201:end,:) = 0;

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_pro_vol_left,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_pro_vol_left,1),size(BW_pro_vol_left,2),size(BW_pro_vol_left,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,3)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.1])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_2.BeamName 'BottomLeft.stl'])

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_pro_vol_top_left,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_pro_vol_top_left,1),size(BW_pro_vol_top_left,2),size(BW_pro_vol_top_left,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,4)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.01])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_2.BeamName 'TopLeft.stl'])

        % tmp code to do only half block
        BW_pro_vol_right = BW_pro_vol;
        BW_pro_vol_top_right = BW_pro_vol_top;
        BW_pro_vol_right(:,1:200,:) = 0;
        BW_pro_vol_top_right(:,1:200,:) = 0;

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_pro_vol_right,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_pro_vol_right,1),size(BW_pro_vol_right,2),size(BW_pro_vol_right,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,5)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.1])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_2.BeamName 'BottomRight.stl'])

        % saving bottom lung blocks as stl file
        [faces,vertices] = extractIsosurface(BW_pro_vol_top_right,0.05);
        triInstrinsic = triangulation(double(faces),double(vertices));
        I = vertices(:,1);
        J = vertices(:,2);
        K = vertices(:,3);
        R = imref3d([size(BW_pro_vol_top_right,1),size(BW_pro_vol_top_right,2),size(BW_pro_vol_top_right,3)]);
        [X,Y,Z] = intrinsicToWorld(R,I,J,K);
        verticesPatientCoords = [X Y Z];

        subplot(3,2,6)     
        trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3)),grid on
        pbaspect([1 1 0.01])

        stlwrite(triInstrinsic,[path H.BeamSequence.Item_2.BeamName 'TopRight.stl'])

        disp('Code Completed');
    end
end
