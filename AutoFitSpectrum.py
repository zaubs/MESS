import imageio, math
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import (r2_score, mean_absolute_error)
from scipy.signal import savgol_filter

image = imageio.imread('TestSpectrum1.png', as_gray=True)

fig, ax = plt.subplots(figsize=(10,6))
# ax = fig.add_subplot(111)


y,x = np.indices(image.shape)



mae = []
scores = []
roll_ransac = []
m_ransac = []
b_ransac = []

for i in range(1,200,5):

	ax.set_xlim(0,500)
	ax.set_ylim(0,500)
	axins = ax.inset_axes([0.0,0.75,0.3,0.25] )#Left, Bottom, Width, Height

	valid_z = (y.ravel()>0) & (image.ravel()>(50+i))
	x_valid = x.ravel()[valid_z]
	y_valid = y.ravel()[valid_z]
	z_valid = image.ravel()[valid_z]
	print(len(z_valid))

	if len(z_valid) > 100:
		sc1 = ax.scatter(x_valid,y_valid, color='yellowgreen', marker='.')
		sc2 = ax.scatter(x_valid,y_valid, color='gold', marker='.')

		ransac = RANSACRegressor(residual_threshold=5).fit(x_valid.reshape(-1,1), y_valid.reshape(-1,1), sample_weight=z_valid**5)
		# ransac.fit(x_valid.reshape(-1,1), y_valid.reshape(-1,1), sample_weight=z_valid**2)
		inlier_mask = ransac.inlier_mask_
		outlier_mask = np.logical_not(inlier_mask)

		line_X = np.arange(x_valid.min(), x_valid.max())[:,np.newaxis]
		line_y_ransac = ransac.predict(line_X)

		prediction = ransac.predict(x_valid.reshape(-1,1))
		# print(ransac.score(x_valid.reshape(-1,1), y_valid.reshape(-1,1)))
		# print(mean_absolute_error(y_valid,prediction))

		mae.append(mean_absolute_error(y_valid,prediction))
		scores.append(ransac.score(x_valid.reshape(-1,1), y_valid.reshape(-1,1)))
		# ax.scatter(x_valid, y_valid, c=z_valid, alpha=0.2, s=20, edgecolor='none', cmap=plt.cm.jet)
		# ax.scatter(x_valid[inlier_mask], y_valid[inlier_mask], color="yellowgreen", marker='.')
		# ax.scatter(x_valid[outlier_mask], y_valid[outlier_mask], color='gold', marker='.')

		sc1.set_offsets(np.c_[x_valid[inlier_mask], y_valid[inlier_mask]])
		sc2.set_offsets(np.c_[x_valid[outlier_mask], y_valid[outlier_mask]])


		z = np.polyfit(x_valid,y_valid, w=z_valid, deg=1)
		p = np.poly1d(z)

		x_plot = np.linspace(x_valid.min(), x_valid.max(), 100)
		y_plot = p(x_plot)

		ax.plot(x_plot, y_plot, '-r', lw=2, label='Poly fit')

		ax.plot(line_X, line_y_ransac, label='RANSAC fit')

		ax.legend(loc='upper right')

		roll_ransac.append(math.degrees(math.atan2((line_y_ransac.max()-line_y_ransac.min()),(line_X.max()-line_X.min()))))
		b_ransac.append(line_y_ransac.max()-((line_y_ransac.max()-line_y_ransac.min())/(line_X.max()-line_X.min()))*line_X.max())
		axins.plot(mae)
	else:
		break
	# print(line_y_ransac.max()-((line_y_ransac.max()-line_y_ransac.min())/(line_X.max()-line_X.min()))*line_X.max())
	plt.axis('equal')

	fig.canvas.draw_idle()
	plt.pause(0.01)
	ax.cla()

	# plt.show()
axins.axvline(10)
plt.waitforbuttonpress()

print(roll_ransac[np.argmin(savgol_filter(mae,21,2))])

plt.figure()
plt.plot(savgol_filter(mae,21,2))
plt.show()