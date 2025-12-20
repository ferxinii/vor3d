clear all

[v, f] = loadHull("build/test_bp.m");

%%

trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;


%%

clear all; close all; hold on;

files = dir("build/*.m");
%files = files(2:end,:);

cmap = lines(numel(files));
for k = 1:numel(files)
    [V, F] = loadHull(fullfile(files(k).folder, files(k).name));
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceAlpha', 0.5, 'FaceColor', cmap(k,:));
    %pause;
end

axis equal; view(3);

%%
function [vertices, faces] = loadHull(filename)
    run(filename);
end
