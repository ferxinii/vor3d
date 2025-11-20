
%% Vcells
[b, fb] = loadHull("build/before.m")
[a, fa] = loadHull("build/after.m")
[bp, fbp] = loadHull("build/bp.m")

trisurf(fb, b(:,1), b(:,2), b(:,3), "FaceColor", "cyan", "FaceAlpha", 0.3)
hold on
trisurf(fa, a(:,1), a(:,2), a(:,3), "FaceColor", "red", "FaceAlpha", 0.3)
trisurf(fbp, bp(:,1), bp(:,2), bp(:,3), "FaceColor", "red", "FaceAlpha", 0.3)

%%
function [vertices, faces] = loadHull(filename)
    run(filename);
end