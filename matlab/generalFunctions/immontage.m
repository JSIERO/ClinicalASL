function displayeddata = immontage(data, rowscolums, noimage)
% ClinicalASL toolbox 2023, JCWSiero
if( nargin == 1 )
    rowscolums = [];
end

displayeddata = MakeMontage(data, rowscolums); % see lines below

if( nargin <= 2 ) % plot the data
    % Plot the image.
    imagesc(displayeddata);

    % Set the colormap.
    colormap(gray(256)); % always use this as default

    % Set the axis.
    axis off;
    axis image;

   set(gcf,'units', 'inches', 'paperposition', [1 1 2.5 2.5]);

end

    function displayeddata = MakeMontage(data, rowscolums)
        %  Check to see if it is something which is shaped to be passed into montage displaty
        if( size(data, 3) == 1 )
            data = squeeze(data);
        end

        if( isempty(rowscolums) )
            %  Calculate the number of columns
            nr = ceil( sqrt(size(data,3)) );
            nc = ceil(size(data,3) / nr);
        else
            nr = rowscolums(1);
            nc = rowscolums(2);
        end

        sx = size(data, 1);
        sy = size(data, 2);

        displayeddata = zeros( sx * nr, sy * nc );

        sl=1;
        for rr = 1:nr
            for cc = 1:nc

                if( sl <= size(data, 3) )
                    row_start = (rr - 1) * sx + 1;
                    row_end = rr * sx;
                    col_start = (cc - 1) * sy + 1;
                    col_end = cc * sy;

                    displayeddata(row_start:row_end, col_start:col_end) = data(:,:,sl);

                    sl = sl + 1;
                end
            end
        end
    end
end