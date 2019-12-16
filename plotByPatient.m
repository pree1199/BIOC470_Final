function plotByPatient(Y,p1se,p2se,p3se,p4se,p5se)
%PLOTBYPATIENT Creates 6 plots of the tSNE visualization with different
%formats of color coding by patient
%color code each point/cell in tSNE by the patient source by using the
%input start and end idices of each patient
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'r','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'b','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'g','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'y','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'k','filled')
legend('Patient 1','Patient 2','Patient 3','Patient 4','Patient 5')
hold off;
title('tSNE Visualization Colored By Source Donor');
xlabel('tSNE_1');
ylabel('tSNE_2');
%color code the tSNE by focusing on individual patients in case the code
%above does not give a clear picture
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'r','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'b','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'b','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'b','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'b','filled')
legend
hold off;
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'b','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'r','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'b','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'b','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'b','filled')
legend
hold off;
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'b','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'b','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'r','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'b','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'b','filled')
legend
hold off;
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'b','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'b','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'b','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'r','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'b','filled')
legend
hold off;
figure;
hold on;
scatter(Y(p1se(1)-1:p1se(2)-1,1),Y(p1se(1)-1:p1se(2)-1,2),5,'b','filled')
scatter(Y(p2se(1)-1:p2se(2)-1,1),Y(p2se(1)-1:p2se(2)-1,2),5,'b','filled')
scatter(Y(p3se(1)-1:p3se(2)-1,1),Y(p3se(1)-1:p3se(2)-1,2),5,'b','filled')
scatter(Y(p4se(1)-1:p4se(2)-1,1),Y(p4se(1)-1:p4se(2)-1,2),5,'b','filled')
scatter(Y(p5se(1)-1:p5se(2)-1,1),Y(p5se(1)-1:p5se(2)-1,2),5,'r','filled')
legend
hold off;
end

