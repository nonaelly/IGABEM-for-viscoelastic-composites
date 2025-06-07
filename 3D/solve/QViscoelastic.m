function [QGi, QKi, Q] = QViscoelastic(QGi, QKi, U1, U2, HG, HK, deltaT, eTauI, i, mu, K, tau, visG, visK)
	QGi = eTauI'.*QGi + (visG.*(1-(tau/deltaT-1).*(1-eTauI)))'.*(HG*U1)...
        -(i>1)*(visG.*(1-tau/deltaT.*(1-eTauI)))'.*(HG*U2);
	QKi = eTauI'.*QKi + (visK.*(1-(tau/deltaT-1).*(1-eTauI)))'.*(HK*U1)...
        -(i>1)*(visK.*(1-tau/deltaT.*(1-eTauI)))'.*(HK*U2);
    Q = 2* mu *sum(QGi,2) + K * sum(QKi,2);
end