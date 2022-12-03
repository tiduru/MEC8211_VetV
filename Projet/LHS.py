import lhsmdu
import matplotlib.pyplot as plt
import numpy as np

lhs_2d = lhsmdu.sample(2,10)
lhs_2d = np.array(lhs_2d)
print("LHS 2D T and rho")
for i in range(len(lhs_2d[0])):
    print(lhs_2d[:,i])

lhs_3d = lhsmdu.sample(3,10)
lhs_3d = np.array(lhs_3d)
print("LHS 3D T rho and b")
for i in range(len(lhs_3d[0])):
    print(lhs_3d[:,i])

# 2D plot
fig1 = plt.figure()
ax1 = fig1.gca()
ax1.set_xticks(np.arange(0,1,0.1))
ax1.set_yticks(np.arange(0,1,0.1))
ax1.scatter(lhs_2d[0], lhs_2d[1])
ax1.set_title(r"Localisation de l'échantillonage" "\n" r"LHS à 2 variables")
ax1.set_xlabel(r"Tension T")
ax1.set_ylabel(r"Densité linéique $\rho$")
ax1.grid()
plt.show()

# 3D plot
fig2 = plt.figure()
ax2 = fig2.add_subplot(projection = "3d")
ax2.set_xticks(np.arange(0,1,0.1))
ax2.set_yticks(np.arange(0,1,0.1))
ax2.scatter(lhs_3d[0], lhs_3d[1], lhs_3d[2])
ax2.set_title(r"Localisation de l'échantillonage" "\n" r"LHS à 3 variables")
ax2.set_xlabel(r"Tension T")
ax2.set_ylabel(r"Densité linéique $\rho$")
ax2.set_zlabel(r"Facteur d'amortissement b")
ax2.grid()
plt.show()