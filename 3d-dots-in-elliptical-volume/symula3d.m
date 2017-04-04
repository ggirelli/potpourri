classdef symula3d
    %SYMULA3D allows to distribute 3D dots in a given volume
    %   This class contains functions that allow to distribute uniformly 2D
    %   dots inside an ellipsoid.
    %   
    %   USE EXAMPLE:
    %   sym = symula3d
    %   sym.N = 1000;
    %   sym.ellSize = [10 7 4];
    %   sym.dr = .25;
    %   sym = sym.run;
    %
    %   Can be coupled with dotterTest using s.symula=sym;
    
    properties
        N = 100;                % Number of points/dots
        ellSize = [5 5 5];      % ellipsoid half-sizes
        vd = [0.13 0.13 0.2];   % voxel size vector [width length height]
        dr = .3;                % dot radius
        imSize = [0 0 0];       % image sizes
        im;                     % image stack
        sd;                     % Gaussian profile SD
        span;                   % camera half-size (3*SD)
        offset = 5;             % image-volume offset [px]s
        uom = 'um';             % units of measure (used in prompts)
        fname = 'test.tif';    % output image file name
        write = 1;              % boolean to write output TIFF
        display = 1;            % boolean to display the result
        points;                 % dot centers in cartesian coordinates
        intensityLim;           % pixel intensity limits
        intensityList;          % pixel intensities
    end
    
    methods
        function instance = symula3d()
            % Remember to change the default properties before hitting "run"!
            display('Remember to change the default properties before hitting "run"!');
            instance.imSize = ceil(2 * instance.ellSize ./ instance.vd);            % image sizes
            instance.imSize = instance.imSize  + 2 * instance.offset;
            instance.im = zeros(instance.imSize);                                   % black image stack
            instance.sd = instance.dr / 2.35482;                                    % Gaussian sigma
            instance.span = ceil(max(4 * instance.sd ./ instance.vd));              % half-size of the camera
        end
        function instance = run(instance)
            % Runs the pipeline with the instance properties.
            
            tic;
            
            disp('> Calculating variables.');
            instance.imSize = ceil(2 * instance.ellSize ./ instance.vd);            % image sizes
            instance.imSize = instance.imSize + 2 * instance.offset;
            instance.im = zeros(instance.imSize);                                   % black image stack
            instance.sd = instance.dr / 2.35482;                                    % Gaussian sigma
            instance.span = ceil(max(4 * instance.sd ./ instance.vd));              % half-size of the camera
            
            disp('> Retrieve points. (the origin is located at the center of the ellipsoid)');
            instance.points = symula3d.uniformPointsInEllipsoid(instance.N,...
                instance.ellSize(1), instance.ellSize(2), instance.ellSize(3));
            
            disp('> Traslate origin to the bottom left corner of the stack and add the specified offset.');
            instance.points(:,1) = instance.points(:,1) + instance.ellSize(1) + instance.offset * instance.vd(1);
            instance.points(:,2) = instance.points(:,2) + instance.ellSize(2) + instance.offset * instance.vd(2);
            instance.points(:,3) = instance.points(:,3) + instance.ellSize(3) + instance.offset * instance.vd(3);
            
            disp('> Calculate pixel intensity.');
            instance.intensityLim = [1 0];
            for n = 1:instance.N
                pixel = symula3d.cartesianToPixel(instance.points(n,:), instance.vd);
                camera = symula3d.getCamera(pixel, instance.span,...
                    'xlim', [1 instance.imSize(1)],...
                    'ylim', [1 instance.imSize(2)],...
                    'zlim', [1 instance.imSize(3)]);
                
                % Pixel-to-Cartesian
                camera(:,1) = camera(:,1) .* instance.vd(1) - instance.vd(1)/2;
                camera(:,2) = camera(:,2) .* instance.vd(2) - instance.vd(2)/2;
                camera(:,3) = camera(:,3) .* instance.vd(3) - instance.vd(3)/2;
                
                intensities = symula3d.getCameraIntensity(instance.points(n,:), camera, 'sd', instance.sd);
                instance.intensityLim(1) = min([instance.intensityLim(1); intensities(:,4)]);
                instance.intensityLim(2) = max([instance.intensityLim(2); intensities(:,4)]);
                
                % Cartesian-to-Pixel
                intensities(:,1) = ceil(intensities(:,1) / instance.vd(1));
                intensities(:,2) = ceil(intensities(:,2) / instance.vd(2));
                intensities(:,3) = ceil(intensities(:,3) / instance.vd(3));
                
                for k = 1:size(intensities,1)
                    cur = intensities(k,:);
                    instance.im(cur(1), cur(2), cur(3)) = instance.im(cur(1), cur(2), cur(3))+cur(4);
                end
                
                instance.intensityList = [instance.intensityList; intensities];
            end
            
            disp('> Normalize instance.intensityList.');
            instance.im = instance.im - instance.intensityLim(1);
            instance.im(instance.im < 0) = 0;
            instance.im = instance.im / instance.intensityLim(2);
            
            if instance.write
                disp('> Write output.');
                imwrite(instance.im(:,:,1), instance.fname, 'Compression', 'none');
                for K = 2:size(instance.im, 3)
                    imwrite(instance.im(:, :, K), instance.fname,...
                        'WriteMode', 'append',  'Compression','none');
                end
            end
            
            if instance.display
                disp('> Visualize result.');
                isosurface(instance.im);
                view(3);
                axis equal;
            end
            
            disp('> Done.');
            
            toc;
        end
    end
    methods (Static)
        function points = uniformPointsInCuboid(N, a, b, c)
            % Uniformly distributes N points in a cuboid
            % of provided dimensions. The origin is located in the center of the cuboid.
            %
            % points = uniformPointsInCuboid(N, a, b, c)
            % N is the number of points.
            % a, b and c are the half-dimensions of the cuboid.
            
            points = zeros(N, 3);
            points(:,1) = rand(N, 1) * 2 * a - a;
            points(:,2) = rand(N, 1) * 2 * b - b;
            points(:,3) = rand(N, 1) * 2 * c - c;
        end
        function points = uniformPointsInEllipsoid(N, a, b, c)
            % Uniformly distribute N points in an ellipsoid
            % of provided dimensions. The origin is located in the center of the ellipsoid.
            %
            % points = uniformPointsInEllipsoid(N, a, b, c)
            % N is the number of points.
            % a, b and c are the half-dimensions of the ellipsoid.
            
            points = nan(N, 3);
            while 0 ~= size(find(isnan(points)))
                candidates = symula3d.uniformPointsInCuboid(N, a, b, c);
                ecd = (candidates(:,1)/a).^2 + (candidates(:,2)/b).^2 + (candidates(:,3)/c).^2;
                
                % Add good candidates to output
                toadd = candidates(ecd <= 1,:);
                ind = find(isnan(points(:,1)), 1);
                if ind+size(toadd, 1) >= N
                    points(ind:end,:) = toadd(1:(N-ind+1),:);
                else
                    points(ind:(ind+size(toadd, 1)-1),:) = toadd;
                end
            end
            points = points(1:N, :);
        end
        function intensity = normFromDist(dc, pc, varargin)
            % Calculates the intensity of a pixel, based on the distance of its center from the dot center.
            % 
            % intensity = normFromDist(dc, pc)
            % dc and pc are, respectively, the dot and pixel center coordinates [x y z].
            %
            % intensity = normFromDist(dc, pc, 'mean', 0)
            % intensity = normFromDist(dc, pc, 'sd', 1)
            % intensity = normFromDist(dc, pc, 'mean', 0, 'sd', 1)
            % sd and mean are optional.
            
            % Read optional parameters
            mean = 0; sd = 1;
            for n = 1:2:numel(varargin)
                if strcmp('mean', varargin{n})
                    mean = varargin{n+1};
                elseif strcmp('sd', varargin{n})
                    sd = varargin{n+1};
                end
            end
            
            % Calculate intensity
            distance = sqrt(sum((dc - pc).^2));
            intensity = 1 / (sd * sqrt(2 * pi)) * exp(-(distance - mean)^2 / (2 * sd^2));
        end
        function grid = getCamera(pc, span, varargin)
            % Provides the coordinates of the pixels around a certain point
            % 
            % grid = getCamera(pc, span)
            % pc [x y z] are the coordinates of the pixel in the center of the square camera grid.
            % span is the half-size of the the square camera (in px).
            %
            % grid = getCamera(pc, span, 'xlim', [0 1])
            % grid = getCamera(pc, span, 'ylim', [0 1])
            % grid = getCamera(pc, span, 'zlim', [0 1])
            % grid = getCamera(pc, span, 'xlim', [0 1], 'ylim', [0 1], 'zlim', [0 1])
            % xlim, ylim and zlim set the lower and upper limit for x, y and z coordinates of camera pixels.
            
            % Build grid (x and y coordinates)
            [X, Y, Z] = meshgrid((pc(1)-span):(pc(1)+span),...
                (pc(2)-span):(pc(2)+span),...
                (pc(3)-span):(pc(3)+span));
            grid = ([reshape(X, 1, numel(X));...
                reshape(Y, 1, numel(Y));...
                reshape(Z, 1, numel(Z))])';
            
            % Read optional parameters
            xlim = NaN; ylim = NaN; zlim = NaN;
            for n = 1:2:numel(varargin)
                if strcmp('xlim', varargin{n})
                    xlim = varargin{n+1};
                elseif strcmp('ylim', varargin{n})
                    ylim = varargin{n+1};
                elseif strcmp('zlim', varargin{n})
                    zlim = varargin{n+1};
                end
            end
            
            % Apply x and y limits
            if ~isnan(xlim)
                grid = grid(grid(:,1) >= xlim(1) & grid(:,1) <= xlim(2), :);
            end
            if ~isnan(ylim)
                grid = grid(grid(:,2) >= ylim(1) & grid(:,2) <= ylim(2), :);
            end
            if ~isnan(zlim)
                grid = grid(grid(:,3) >= zlim(1) & grid(:,3) <= zlim(2), :);
            end
        end
        function intensities = getCameraIntensity(dc, grid, varargin)
            % Calls normFromDist for every pixel in grid.
            %
            % intensities = getCameraIntensity(dc, grid)
            % dc contains the dot center coordinates [x y z].
            % grid is an nx3 matrix that contains the coordinates of the grid pixel center.
            %
            % intensities = getCameraIntensity(dc, grid, 'mean', 0)
            % intensities = getCameraIntensity(dc, grid, 'sd', 1)
            % intensities = getCameraIntensity(dc, grid, 'mean', 0, 'sd', 1)
            % sd and mean are optional.
            
            % Read optional parameters
            mean = 0; sd = 1;
            for n = 1:2:numel(varargin)
                if strcmp('mean', varargin{n})
                    mean = varargin{n+1};
                elseif strcmp('sd', varargin{n})
                    sd = varargin{n+1};
                end
            end
            
            % Calculate intensities
            intensities = zeros(size(grid,1), 4);
            for n = 1:size(grid, 1)
                intensities(n,:) = [grid(n,:) symula3d.normFromDist(dc, grid(n,:),...
                    'mean', mean, 'sd', sd)];
            end
            
            % Normalize intensities
            intensities(:, 4) = intensities(:, 4) - min(intensities(:, 4));
            intensities(intensities(:, 4) < 0, 4) = 0;
            intensities(:, 4) = intensities(:, 4) / max(intensities(:, 4));
        end
        function coords = cartesianToPixel(dc, vd)
            % Converts cartesian coordinates to pixel coordinates.
            %
            % coords = cartesianToPixel(dc, vd)
            % dc are the cartesian coordinates of the dot center.
            % vd are the voxel dimensions.
            
            coords = ceil(dc ./ vd);
        end
        function coords = pixelToCartesian(dc, vd)
            % Converts pixel coordinates to cartesian coordinates
            %
            % coords = pixelToCartesian(dc, vd)
            % dc are the cartesian coordinates of the dot center.
            % vd are the voxel dimensions.
            
            coords = dc .* vd - vd/2;
        end
    end
    
end

