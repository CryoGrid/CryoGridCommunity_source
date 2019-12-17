% figures new

% 3
i = 1;
for i = 1:368 %360 for 4
plot(out.STRATIGRAPHY{1, i}{1, 1}.STATVAR.T(1:end))
hold on
end

% %% 4
% i = 1;
% for i = 1:360 
% plot(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.T(1:end))
% plot(out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(1:end))
% hold on
% end
% 
% %% 1
% i = 1;
% for i = 1:360 
% plot(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.T(1:end))
% plot(out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(1:end))
% plot(out.STRATIGRAPHY{1, i}{3, 1}.STATVAR.T(1))
% hold on
% end

% %% 2
% i = 1;
% for i = 1:368
% plot(out.STRATIGRAPHY{1, 1}{1, 1}.STATVAR.T(1:end))
% plot(out.STRATIGRAPHY{1, i}{2, 1}.STATVAR.T(1:end))
% hold on
% end