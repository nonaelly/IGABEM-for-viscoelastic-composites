function [QGi, QKi, Q] = QViscoelastic(QGi, QKi, U1, U2, HG, HK, deltaT, eTaoI, i, mu, K, tau, visG, visK)
	QGi = eTaoI'.*QGi + (visG.*(1-(tau/deltaT-1).*(1-eTaoI)))'.*(HG*U1)...
        -(i>1)*(visG.*(1-tau/deltaT.*(1-eTaoI)))'.*(HG*U2);
	QKi = eTaoI'.*QKi + (visK.*(1-(tau/deltaT-1).*(1-eTaoI)))'.*(HK*U1)...
        -(i>1)*(visK.*(1-tau/deltaT.*(1-eTaoI)))'.*(HK*U2);
    Q = 2* mu *sum(QGi,2) + K * sum(QKi,2);
end