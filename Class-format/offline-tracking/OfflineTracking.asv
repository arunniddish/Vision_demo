classdef OfflineTracking

    properties

        number_of_markers;
        centroids;
        PrevPt;
        P0;
        CurrPt;
        cent;
        
        vread;
        numberOfFrames;
        vwrite;

    end

    methods

        function obj = OfflineTracking(params)

            obj.number_of_markers = params.number_of_markers;
            obj.PrevPt = [];
            obj.P0 = [];
            obj.CurrPt = [];
            obj.cent = [];

            obj.vread = params.vread;
            obj.numberOfFrames = obj.vread.NumberOfFrame;
            obj.vwrite = params.vwrite;

            obj.centroids = zeros(obj.numberOfFrames,3*obj.number_of_markers);

        end

        function tracking_data = tracking(obj)

            for k = 1:obj.numberOfFrames

                thisFrame = read(obj.vread,k);
                newim = createMaskCarpetBlue(thisFrame);
                newim = bwareaopen(newim,20);
                newim = imfill(newim, 'holes');
%                 axis on;  % To be changed
                [labeledImage, numberOfRegions] = bwlabel(newim);
                count = 0;
%                 cent = [];
                obj.cent = zeros(numberOfRegions,2);

                stats = regionprops(labeledImage, 'BoundingBox','Centroid','Area','EquivDiameter');

                 for rb = 1:numberOfRegions
                     count = count + 1;
                     obj.cent(count,:) = stats(rb).Centroid;
                     obj.cent(count,2) = 1080 - obj.cent(count,2) ;  % Correction for y-axis.
                 end

                  zc = zeros(size(obj.cent,1),1);
                  obj.cent = [obj.cent,zc];

                  if k == 1
                    obj.P0 = obj.cent;
                    obj.PrevPt = obj.cent;
                    obj.centroids = data_logging(obj,k);
                  end

                  if k ~= 1

                     obj.centroids = nearest_neighbor(obj,count);
                     obj.centroids = data_logging(obj,count);
%                      [~,~,theta(k,:),trans(k,:)] = pose_estimation(obj.CurrPt,obj.PrevPt,k);
                     [Rot,T] = pose_estimation(obj,obj.centroids,obj.PrevPt,k);
%                      theta(k,:) = reshape(Rot,[1,9]);
%                      trans(k,:) = T';
%                      [~,~,theta_G(k,:),trans_G(k,:)] = pose_estimation(obj.P0,obj.PrevPt,k);
                     [Rot,T] = pose_estimation(obj.P0,obj.PrevPt,k);
%                      theta_G(k,:) = reshape(Rot,[1,9]);
%                      trans_G(k,:) = T';
                     obj.PrevPt = obj.CurrPt;
                  end
                  
            end
            tracking_data = cat(2,obj.centroids,theta,trans,theta_G,trans_G);
        end

        function centroids = nearest_neighbor(obj,count)      %Change the name of the centroid ->
            
            obj.CurrPt = zeros(obj.number_of_markers,3);

            if(count > obj.number_of_markers)
                for i = 1:obj.number_of_markers
                    for j = 1:count
                        X = [obj.PrevPt(i,:);obj.cent(j,:)];
                        d(j) = pdist(X,'euclidean');
                    end
                    [dmin,ind] = min(d);  
                    if(dmin < 15)
                        obj.CurrPt(i,:) = obj.cent(ind,:);
                    end
                end
            end

            if(count <= obj.number_of_markers)
                for i = 1:count
                    for j = 1:obj.number_of_markers
                        X = [obj.cent(i,:);obj.PrevPt(j,:)];
                        d(j) = pdist(X,'euclidean');
                    end
                    [dmin,ind] = min(d);
                    if(dmin < 15)
                        obj.CurrPt(ind,:) = obj.cent(i,:);
                    end
                end
            end
            
            clear d;
            TF = obj.CurrPt(:,1);  % Writing the 1st column of resrvd
            index = find(TF == 0);  % Finding those rows which is empty
            val = isempty(index);   % Checking whether the index is empty  
            
            if(val == 0)
                centroids = occlusion(obj,index);   %Change the name of the centroid ->
            end
            if(val~=0)
                centroids = obj.CurrPt;       %Change the name of the centroid ->
            end
        end 

        function centroids = occlusion(obj,index)   %Change the name of the centroid ->

            newPrevPt = obj.PrevPt;
            newP0 = obj.P0; 
            newPrevPt(index(1:size(index,1)),:) = 0;
            newP0(index(1:size(index,1)),:) = 0;
%             [Rot,T,~,~] = pose_estimation(newPrevPt,obj.CurrPt); % SE2 w.r.t previous frame
            [Rot,T] = pose_estimation(obj,newPrevPt,obj.CurrPt); % SE2 w.r.t previous frame

            for gg = 1:size(index,1)
                newPt = Rot*(obj.PrevPt(index(gg),:))' + T;
                obj.CurrPt(index(gg),:) = newPt;
            end

            centroids = obj.CurrPt;

        end

        function centroids = data_logging(obj,k)
            for i = 1:obj.number_of_markers
                centroids(k,(3*i)-2:(3*i)) = obj.PrevPt(i,:);
            end
        end

        function [Rot,T,theta(k,:),trans(k,:)] = pose_estimation(obj,A,B,k)
%         function [Rot,T] = pose_estimation(obj,A,B,k)
            [Rot,T] = rigid_transform_3D(A',B');  % SE2 w.r.t previous frame
            theta(k,:) = reshape(Rot,[1,9]);   % Changed ANM
            trans(k,:) = T';
        end

        function plot(obj)
            figure(10)
            imshow(thisFrame)
            set(gcf, 'Position',  [100, 100, 1000, 1000])
            hold on
            plot(obj.PrevPt(:,1),obj.PrevPt(:,2),'g*','LineWidth',0.5,'MarkerSize',2)
            caption = sprintf('%d blobs found in frame #%d 0f %d', count, k, obj.numberOfFrames);
            title(caption, 'FontSize', 20);
            axis on;
            hold off
            pframe = getframe(gcf);
            writeVideo(obj.vwrite,pframe);
        end
    
    end
end