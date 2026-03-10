clear all


%%
P = readmatrix("problematic_sphere.csv");

scatter3(P(:,1), P(:,2), P(:,3), "green"); hold on;
P(:,4) = P(:,4)+1;

%%
[v, f] = loadHull("build/hole.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

P = readmatrix("build/nonmani.txt");
scatter3(v(P+1,1), v(P+1,2), v(P+1,3), "green"); hold on;



%%

[v, f] = loadHull("build/prev.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

P = readmatrix("build/nonmani.txt");
scatter3(v(P+1,1), v(P+1,2), v(P+1,3), "green"); hold on;

%%
[v, f] = loadHull("build/deleted.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

%%

k = readmatrix("build/visible.txt");
scatter3(v(k(:)+1,1), v(k(:)+1,2), v(k(:)+1,3), 'blue'); hold on;

%%

h = readmatrix("build/horizon.txt");
scatter3(v(h(:)+1,1), v(h(:)+1,2), v(h(:)+1,3), 'yellow'); hold on;

%%
a = 197;
scatter3(v(a,1), v(a,2), v(a,3), 'red'); 

%%
P = readmatrix("build/test_circum.txt");

scatter3(P(2:end,1), P(2:end,2), P(2:end,3)); hold on;
scatter3(P(1,1), P(1,2), P(1,3), "red");

%%
[vb, fb] = loadHull("build/before.m");
[va, fa] = loadHull("build/after.m");

trisurf(fb, vb(:,1), vb(:,2), vb(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;
trisurf(fa, va(:,1), va(:,2), va(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);


P = readmatrix("build/test_intersect.txt");
for i = 1:size(P, 1)
    scatter3(P(i,1), P(i,2), P(i,3));
end

%%

clear all; close all; %hold on;

files = dir("build/v*.m");
cmap = lines(numel(files));
tot_vol = 0
for k = 1:numel(files)
    if (files(k).bytes == 0)
        k
        continue;
    end
    [V, F] = loadHull(fullfile(files(k).folder, files(k).name));
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceAlpha', 0.5, 'FaceColor', cmap(k,:));

    [ch, vol] = convhull(V);
    vol
    tol_vol = tot_vol + vol;

    axis equal; view(3);
    pause;


end

view(3);


a = 5;
h = readmatrix("seeds.csv");
scatter3(h(a,1), h(a,2), h(a,3), 'green'); 



%%
k = 3;
[v, f] = loadHull(fullfile(files(k).folder, files(k).name));
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'Cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

%%
[v, f] = loadHull("v2.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.3, 'LineWidth', 0.001);
hold on;

%%
[v, f] = loadHull("build/v2.m");
trisurf(f, v(:,1), v(:,2), v(:,3), 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'LineWidth', 0.001);





%%
function [vertices, faces] = loadHull(filename)
    run(filename);
end
