%% Parameters 

lambda = 3.8e-3;      % wavelength (m)
w0     = 40e-3;       % beam waist (m)
k      = 2*pi/lambda; % wave number

z0 = pi * w0^2 / lambda;   % Rayleigh distance

% Required propagation distances:
z_list_before = [0, 0.5, 1];        % before reflector
z_list_after  = [1,4,6];            % after reflector (from mirror)
z_list_reflector=[3,4,5];           % mirror locations
f = -4*z0;                           % parabolic reflector focus
c=['r','b','g','k','m','c','g'];    % colors for plotting
axis_limit=300;
%% Simulation Grid 

N = 1024;               % points in x-y
L = 3.0;                % physical grid size (meters)
dx = L/N;               % sampling step
center = N/2 + 1;

z_max_est = N * dx^2 / lambda;

min_dx = (sqrt(2)*pi)/k;
if dx > min_dx
    fprintf('Sampling condition satisfied: dx = %e  >  %e\n', dx, min_dx);
else
    fprintf('Sampling condition may be violated: dx = %e  <=  %e\n', dx, min_dx);
end

x = (-N/2:N/2-1)*dx;
y = x;
[X,Y] = meshgrid(x,y);

% Spatial frequencies
fx = (-N/2:N/2-1)/(N*dx);
fy = fx;
[Fx, Fy] = meshgrid(fx, fy);

kx = 2*pi*Fx;
ky = 2*pi*Fy;

%% Initial Gaussian Field 

U0 = exp(-(X.^2 + Y.^2)/(w0^2));   
x_mm = x*1e3;
idx = find(x_mm >= -200 & x_mm <= 200);

figure(1);
set(gcf,'Position',[100 100 1400 450]);   
sgtitle('Gaussian wave beam at z=0','FontWeight','bold');
    subplot(1,3,1);
        imagesc(x_mm(idx), x_mm(idx), abs(U0(idx,idx)).^2);
        title('Gaussian Beam Intensity at z = 0');
        xlabel('x (mm)');
        ylabel('y (mm)');
        axis image;           
        colorbar;
        colormap('hot');

    subplot(1,3,2);
        plot(x*1e3, abs(U0(center,:)).^2,'r','LineWidth',2);
        title('Central Slice at y = 0');
        xlabel('x (mm)');
        ylabel('Intensity');
        xlim([-150 150]);
        grid on;
        pbaspect([1.3 1 1]); % adjusts width-to-height ratio

    subplot(1,3,3);
        surf(x_mm(idx), x_mm(idx), abs(U0(idx,idx)).^2);
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('Intensity');
        title('3D Intensity Profile at z = 0');
        shading interp;         % prevent shading for intensity visualization
        axis image;
        colormap hot;
        daspect([1 1.2 0.0014]); % scaling in X,Y,Z axes

%% Frequency-domain representation and Propagation Factor

U0_f = fftshift(fft2(U0));

prop = @(z) exp(-1j*z*(k-((kx.^2 + ky.^2)/(2*k))));

%% Propagation Before Reflection 
figBefore2D = figure;
tL = tiledlayout(1, length(z_list_before), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tL, 'Beam Intensity before Mirror at Different Propagation Distances', 'FontWeight','bold');

Power_before = zeros(1,length(z_list_before));  % initialize power array

for i = 1:length(z_list_before)
    z = z_list_before(i)*z0;
    Uz_f = U0_f .* prop(z);
    Uz   = ifft2(ifftshift(Uz_f));
    
    % Calculate total power before mirror
    Power_before(i) = sum(sum(abs(Uz).^2)) * dx^2;

    figure(figBefore2D);
    ax = nexttile(i);

    % Plot intensity
    imagesc(ax, x_mm(idx), x_mm(idx), abs(Uz(idx,idx)).^2);
    axis(ax, 'image');         
    xlabel(ax, 'x (mm)');
    ylabel(ax, 'y (mm)');
    title(ax, ['At z = ', num2str(z_list_before(i), '%.1f'), ' z_0']);
    colormap(ax, 'hot');
    colorbar(ax);



    if(i==1)
        Q = length(findobj('Type','figure')); % count existing figures
        figure(Q + 1);                        % open the next figure
    else
        figure(Q+1);
    end
        axis tight;
        plot(x*1e3, abs(Uz(center,:)).^2,c(i),'LineWidth',2);hold on;
        title('Slice of Relflectled Wave Intensity at at y=0 & different multiples of z_o');
        xlabel('x (mm)'); ylabel('Intensity');
        xlim([-150 150]);      
        grid on;
end
legend('z=0','z=0.5z_o','z=z_o');
%% Mirror Reflection and propagation

M = exp(-1j * k * (X.^2 + Y.^2) / (2*f));
% preallocations for the figure and power arrays
Power_after = zeros(length(z_list_reflector), length(z_list_after));
figAfter2D = gobjects(1, length(z_list_reflector));
%Uout_f_storage=zeros(length(z_list_reflector));
for j=1:length(z_list_reflector)
    z = z_list_reflector(j) * z0; 
    Uin_f = U0_f .* prop(z);
    Uin   = ifft2(ifftshift(Uin_f));
    Uout = M .* Uin;
    Uout_f = fftshift(fft2(Uout));
    %Uout_f_storage(j) = Uout_f; % Store the reflected fiel

    figAfter2D(j) = figure;
    tL_after = tiledlayout(1, length(z_list_after), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tL_after, ['Reflected Beam Intensities after Mirror at z_{mirror} = ', num2str(z_list_reflector(j)), ' z_0'], ...
        'FontWeight', 'bold');
    for i = 1:length(z_list_after)
        z = z_list_after(i)*z0;
        Uz_f = Uout_f .* prop(z);
        Uz   = ifft2(ifftshift(Uz_f));

        % Calculate total power
        Power_after(j,i) = sum(sum(abs(Uz).^2)) * dx^2;
        
        % subplot(1,length(z_list_after),i);
        figure(figAfter2D(j));
        ax = nexttile(i);
        imagesc(ax, x_mm(idx), x_mm(idx), abs(Uz(idx,idx)).^2);
        axis(ax, 'image');
        xlabel(ax, 'x (mm)');
        ylabel(ax, 'y (mm)');
        title(ax, ['At z = ', num2str(z_list_after(i), '%.1f'), ' z_0']);
        colormap(ax, 'hot');
        colorbar(ax);
    
        if(i==1)
            Q = length(findobj('Type','figure'));   % count existing figures
            fig(j)=figure(Q + 1);                   % open the next figure
            sgtitle(['Slice of Reflected Wave intensity from z_{mirror}=', ...
            num2str(z_list_reflector(j)),' z_o at different distances'], ...
           'FontSize', 12, 'FontWeight', 'bold');
        else
            figure(Q+1);
        end 

        plot(x*1e3, abs(Uz(center,:)).^2,c(i),'LineWidth',2); hold on;
        xlabel('x (mm)'); ylabel('Intensity');
        xlim([-axis_limit axis_limit]);    
        grid on;
    end
legend('At z=z_o','z=4z_o','z=6z_o','FontWeight','bold');
end
%% Slice of Reflected Wave intensity from different mirror positions at different distances from the mirror

figAll_slices = figure;
sgtitle(['Slice of Reflected Wave intensity from different mirror positions at different distances ' ...
    'from the mirror'],'FontWeight','bold');
subplot(3,1,1);   % 3 rows, 1 column, first subplot
copyobj(allchild(fig(1).CurrentAxes), gca);   % copy contents
title('Slice of Reflected Wave intensity from mirror at z=3z_o');
xlim([-axis_limit axis_limit]);


subplot(3,1,2);
copyobj(allchild(fig(2).CurrentAxes), gca);
title('Slice of Reflected Wave intensity from mirror at Z=4z_o');
xlim([-axis_limit axis_limit]);


subplot(3,1,3);
copyobj(allchild(fig(3).CurrentAxes), gca);
title('Slice of Reflected Wave intensity from mirror at z=5z_o');
xlim([-axis_limit axis_limit]);

%% Reflected Beam intensity for different mirror position at different Distances

figuu = figure;
main = tiledlayout(figuu, 3, 1);   % 3 rows, 1 column
title(main, ['Reflected Beam intensity for different mirror position at different' ...
    ' propagating Distances'], 'FontWeight', 'bold');

rowTitles = ["Mirror position at z=3z_o" , "Mirror position at z=4z_o" , "Mirror position at z=5z_o"];

for r = 1:3

    % Create a nested layout for each row (1 row × 3 columns)
    row = tiledlayout(main, 1, 3);
    row.Layout.Tile = r;    % Attach this row to the r-th row of MAIN
    title(row, rowTitles(r), 'FontSize', 14, 'FontWeight', 'bold');

    % Get the 3 axes stored inside your original figures for mirror r
    axSource = flip(findobj(figAfter2D(r), 'Type', 'Axes'));

    for c = 1:3
        ax = nexttile(row);

        % Copy objects from the source subplot into the new tile
        copyobj(allchild(axSource(c)), ax);

        xlabel(ax, 'x (mm)');
        ylabel(ax, 'y (mm)');
        axis(ax, 'image');
        colormap(ax, hot);
        colorbar(ax);

        % Column title
        title(ax, ['Propagation at z = ' num2str(z_list_after(c)) ' z_0']);
    end
end


disp('Power before mirror:');
disp(Power_before);

disp('Power after mirror:');
disp(Power_after); 
%% % --- compute w(z) by locating x where Intensity == 1/e^2 * Imax ---
z_norm_axis = linspace(-3,3,1001);  
w = NaN(size(z_norm_axis));         
Uprop_f_base = U0_f;       

% Correct indices for center (N/2 + 1)
center_col = N/2 + 1;
center_row = N/2 + 1; % Adjusted to be consistent

for s = 1:length(z_norm_axis)
    z_current = z_norm_axis(s)*z0; % Physical propagation distance
    Uz_f = Uprop_f_base .* prop(z_current);
    Uz   = ifft2(ifftshift(Uz_f));

    int_slice = abs(Uz(center_row, :)).^2; 

    I_max = max(int_slice);
    thr   = I_max * exp(-2); % 1/e^2 threshold

    % Find threshold crossing to the right of the peak
    % Note: center_col is now correct (index of x=0)
    right_part = int_slice(center_col:end); 
    j = find(right_part <= thr, 1, 'first');

    if ~isempty(j)
        % Indices for interpolation (must be 1-based, global indices)
        i2 = center_col + j - 1; 
        i1 = i2 - 1;

        % Linear interpolation to find the exact crossing point
        v1 = int_slice(i1);
        v2 = int_slice(i2);
        x1 = x(i1);
        x2 = x(i2);
        t = (thr - v1)/(v2 - v1);
        x_cross = x1 + t*(x2 - x1);

        % Beam width w(s) is the radius (distance from center x=0)
        w(s) = abs(x_cross); 
    end
end

% Theoretical Calculation (w_t)
% Since z_norm_axis is z/z0, the correct formula is used directly.
wt = w0 * sqrt(1 + z_norm_axis.^2); 

%% Plot
figure;

plot(z_norm_axis,  w*1e3, 'b', 'LineWidth', 3); hold on;
plot(z_norm_axis,  wt*1e3, 'r--', 'LineWidth', 2);
plot(z_norm_axis, -w*1e3, 'b', 'LineWidth', 3); hold on;  % mirror (negative side)
plot(z_norm_axis, -wt*1e3, 'r--', 'LineWidth', 2); % mirror theoretical curve

xlabel('Normalized Distance (z / z_0)');
ylabel('Beam radius w (mm)');
title('Beam Radius');
legend('w(z) Sim', 'w(z) Theory');
grid on;

% Final check and output

[~, idx_rayleigh] = min(abs(z_norm_axis - 1));

w_at_z0_sim = abs(w(idx_rayleigh)) * 1e3; 

w_at_z0_theory = w0 * sqrt(2) * 1e3; 


w_min_sim = min(abs(w))*1e3;
w0_theory = w0*1e3;


fprintf('\n--- Beam Waist Results ---\n');
fprintf('Theoretical Waist (w0): %.1f mm\n', w0_theory);
fprintf('Simulated Minimum Waist: %.1f mm\n', w_min_sim);
fprintf('Difference (w0): %.2e mm\n', w_min_sim - w0_theory);

fprintf('\n--- Rayleigh Range (z/z0 = 1) Results ---\n');
fprintf('Theoretical Width (w(z0)): %.1f mm\n', w_at_z0_theory);
fprintf('Simulated Width (w(z0)): %.1f mm\n', w_at_z0_sim);
fprintf('Difference (z0): %.2e mm\n', w_at_z0_sim - w_at_z0_theory);
% z_prop_norm_max = 10; 
% z_prop_norm_axis = linspace(0, z_prop_norm_max, 3001); 
% 
% % Plotting setup
% figure;
% sgtitle('Beam Width Envelope After Parabolic Reflector', 'FontWeight', 'bold');
% legend_entries = {};
% 
% % Correct indices for center (N/2 + 1)
% center_col = N/2 + 1;
% center_row = N/2 + 1; 
% 
% for j = 1:length(z_list_reflector)
% 
%     z_mirror_norm = z_list_reflector(j);
%     Uprop_f_base = Uout_f_storage{j}; % Use the stored reflected field (Uout_f)
% 
%     w_reflected = NaN(size(z_prop_norm_axis));
% 
%     % Loop through propagation distance after the mirror
%     for s = 1:length(z_prop_norm_axis)
%         z_current = z_prop_norm_axis(s) * z0; % Physical propagation distance from mirror
% 
%         Uz_f = Uprop_f_base .* prop(z_current);
%         Uz   = ifft2(ifftshift(Uz_f));
% 
%         int_slice = abs(Uz(center_row, :)).^2; 
% 
%         I_max = max(int_slice);
%         thr   = I_max * exp(-2); % 1/e^2 threshold
% 
%         % Find threshold crossing
%         right_part = int_slice(center_col:end); 
%         k_idx = find(right_part <= thr, 1, 'first');
% 
%         if ~isempty(k_idx)
%             % Interpolation logic
%             i2 = center_col + k_idx - 1; 
%             i1 = i2 - 1;
%             v1 = int_slice(i1);
%             v2 = int_slice(i2);
%             x1 = x(i1);
%             x2 = x(i2);
%             t = (thr - v1)/(v2 - v1);
%             x_cross = x1 + t*(x2 - x1);
%             w_reflected(s) = abs(x_cross);
%         end
%     end
% 
%     % Shift the Z axis for plotting: Z_Plot = Z_Mirror + Z_Prop_After
%     z_plot_axis = z_mirror_norm + z_prop_norm_axis;
%     w_mm = w_reflected * 1e3;
% 
%     % Plotting (using the color 'c' array)
% 
%     % Plot Positive Envelope
%     plot(z_plot_axis, w_mm, 'Color', c(j), 'LineWidth', 2); hold on;
%     % Plot Negative Envelope
%     plot(z_plot_axis, -w_mm, 'Color', c(j), 'LineWidth', 2, 'HandleVisibility', 'off'); 
% 
%     legend_entries{j} = ['Sim. (z_{mir} = ' num2str(z_mirror_norm) ' z_0)'];
% end
% 
% % Plot the propagation axis (Z=0) line
% plot(linspace(0, z_list_reflector(end) + z_prop_norm_max, 100), zeros(1, 100), 'k:', 'HandleVisibility', 'off');
% 
% % Plot the mirror locations
% for j = 1:length(z_list_reflector)
%     z_mirror_norm = z_list_reflector(j);
%     plot([z_mirror_norm, z_mirror_norm], ylim, [c(j) '--'], 'HandleVisibility', 'off'); % Vertical dashed line
% end
% 
% xlabel('Absolute Normalized Distance (z / z_0)');
% ylabel('Beam radius w (mm)');
% legend(legend_entries, 'Location', 'northeast');
% grid on;
% hold off;