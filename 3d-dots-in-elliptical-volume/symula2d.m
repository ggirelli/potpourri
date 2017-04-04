classdef symula2d
    %SYMULA2D allows to distribute dots in a given volume
    %   This class contains functions that allow to distribute uniformly 2D
    %   dots inside an ellipsoid.
    
    properties
        N = 100;                % Number of points/dots
        a = 10;                  % x ellipsoid half-size
        b = 5;                  % y ellipsoid half-size
        c = 1.5;                  % z ellipsoid half-size
        vx = .13;               % x voxel size (width)
        vy = .13;               % y voxel size (length)
        vz = .2;                % z voxel size (height)
        dr = .3;                % dot radius
        vd;                     % voxel size vector
        ix;                     % image x size
        iy;                     % image y size
        iz;                     % image z size (number of slices)
        im;                     % image stack
        sd;                     % Gaussian profile SD
        span;                   % camera half-size (3*SD)
        uom = 'um';             % units of measure (used in prompts)
        fname = 'test.tiff';    % output image file name
        display = 1;         % boolean to display the result
        points;                 % dot centers in cartesian coordinates
        intensityList;          % pixel intensities
    end
    
    methods
        function instance = symula2d()
            % Remember to change the default properties before hitting "run"!
            display('Remember to change the default properties before hitting "run"!');
            instance.vd = [instance.vx instance.vy instance.vz];    % voxel dimensions
            instance.ix = ceil(2 * instance.a / instance.vd(1));        % x image size in vx
            instance.iy = ceil(2 * instance.b / instance.vd(2));        % y image size in vx
            instance.iz = ceil(2 * instance.c / instance.vd(3));        % z image size in vx (i.e.: number of slides)
            instance.im = zeros(instance.ix, instance.iy, instance.iz);                 % black image stack
            instance.sd = instance.dr / 2.35482;               % Gaussian sigma
            instance.span = ceil(max(3 * instance.sd / instance.vd(1), 3 * instance.sd / instance.vd(2)));   % half-size of the camera
        end
        function instance = run(instance)
            % Runs the pipeline with the instance properties.
            
            tic;
            
            disp('> Calculating variables.');
            instance.vd = [instance.vx instance.vy instance.vz];    % voxel dimensions
            instance.ix = ceil(2 * instance.a / instance.vd(1));        % x image size in vx
            instance.iy = ceil(2 * instance.b / instance.vd(2));        % y image size in vx
            instance.iz = ceil(2 * instance.c / instance.vd(3));        % z image size in vx (i.e.: number of slides)
            instance.im = zeros(instance.ix, instance.iy, instance.iz);                 % black image stack
            instance.sd = instance.dr / 2.35482;               % Gaussian sigma
            instance.span = ceil(max(3 * instance.sd / instance.vd(1), 3 * instance.sd / instance.vd(1)));   % half-size of the camera
            
            disp('> Retrieve instance.points. (the origin is located at the center of the ellipsoid)');
            instance.points = symula2d.uniformPointsInEllipsoid(instance.N, instance.a, instance.b, instance.c);
            
            disp('> Traslate origin to the bottom left corner of the stack.');
            instance.points(:,1) = instance.points(:,1) + instance.a;
            instance.points(:,2) = instance.points(:,2) + instance.b;
            instance.points(:,3) = instance.points(:,3) + instance.c;
            
            disp('> Calculate pixel intensity.');
            intensityLim = [1 0];
            for n = 1:instance.N
                pixel = symula2d.cartesianToPixel(instance.points(n,:), instance.vd);
                camera = symula2d.getCamera(pixel(1), pixel(2), instance.span, 'xlim', [1 instance.ix], 'ylim', [1 instance.iy]);
                
                % Pixel-to-Cartesian
                camera(:,1) = camera(:,1) .* instance.vd(1) - instance.vd(1)/2;
                camera(:,2) = camera(:,2) .* instance.vd(2) - instance.vd(2)/2;
                
                intensities = symula2d.getCameraIntensity(instance.points(n,:), camera, 'sd', instance.sd);
                intensities = ([intensities(:,1)'; intensities(:,2)'; repmat(pixel(3), 1, size(intensities, 1)); intensities(:,3)'])';
                intensityLim(1) = min([intensityLim(1); intensities(:,3)]);
                intensityLim(2) = max([intensityLim(2); intensities(:,3)]);
                
                % Cartesian-to-Pixel
                intensities(:,1) = ceil(intensities(:,1) / instance.vd(1));
                intensities(:,2) = ceil(intensities(:,2) / instance.vd(2));
                
                for k = 1:size(intensities,1)
                    cur = intensities(k,:);
                    instance.im(cur(1), cur(2), cur(3)) = instance.im(cur(1), cur(2), cur(3))+cur(4);
                end
                
                instance.intensityList = [instance.intensityList; intensities];
            end
            
            disp('> Normalize instance.intensityList.');
            instance.im = instance.im - intensityLim(1);
            instance.im(instance.im < 0) = 0;
            instance.im = instance.im / intensityLim(2);
            
            disp('> Writing output.');
            imwrite(instance.im(:,:,1), instance.fname, 'Compression', 'none');
            for K = 2:size(instance.im, 3)
                imwrite(instance.im(:, :, K), instance.fname, 'WriteMode', 'append',  'Compression','none');
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
                candidates = symula2d.uniformPointsInCuboid(N, a, b, c);
                ecd = (candidates(:,1)/a).^2 + (candidates(:,2)/b).^2 + (candidates(:,3)/c).^2;
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
            % dc and pc are, respectively, the dot and pixel center coordinates [x y].
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
        function grid = getCamera(x, y, span, varargin)
            % Provides the coordinates of the pixels around a certain point
            % 
            % grid = getCamera(x, y, span)
            % x and y are the coordinates of the pixel in the center of the square camera grid.
            % span is the half-size of the the square camera (in px).
            %
            % grid = getCamera(x, y, span, 'xlim', [0 1])
            % grid = getCamera(x, y, span, 'ylim', [0 1])
            % grid = getCamera(x, y, span, 'xlim', [0 1], 'ylim', [0 1])
            % xlim and ylim set the lower and upper limit for x and y coordinates of camera pixels.
            
            % Build grid (x and y coordinates)
            [X, Y] = meshgrid((x-span):(x+span), (y-span):(y+span));
            grid = ([reshape(X, 1, numel(X));reshape(Y, 1, numel(Y))])';
            
            % Read optional parameters
            xlim = NaN; ylim = NaN;
            for n = 1:2:numel(varargin)
                if strcmp('xlim', varargin{n})
                    xlim = varargin{n+1};
                elseif strcmp('ylim', varargin{n})
                    ylim = varargin{n+1};
                end
            end
            
            % Apply x and y limits
            if ~isnan(xlim)
                grid = grid(grid(:,1) >= xlim(1) & grid(:,1) <= xlim(2), :);
            end
            if ~isnan(ylim)
                grid = grid(grid(:,2) >= ylim(1) & grid(:,2) <= ylim(2), :);
            end
        end
        function intensities = getCameraIntensity(dc, grid, varargin)
            % Calls normFromDist for every pixel in grid.
            %
            % intensities = getCameraIntensity(dc, grid)
            % dc contains the dot center coordinates [x y].
            % grid is an nx2 matrix that contains the coordinates of the grid pixel center.
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
            intensities = zeros(size(grid,1), 3);
            for n = 1:size(grid, 1)
                intensities(n,:) = [grid(n,:) symula2d.normFromDist(dc(1:2), grid(n,:), 'mean', mean, 'sd', sd)];
            end
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
            coords = dc .* vd;
        end
    end
    
end

