import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class DensityOfStatesCalculator:

    def __init__(self, dispersion_relation: callable, grid_points):
        self.dispersion = dispersion_relation
        self.grid = grid_points

    def get_2d_contourf_plot(self, no_of_levels = 50, cmap = 'hot'):
        X, Y = np.meshgrid(self.grid[0], self.grid[1])
        Z = self.dispersion((X, Y))

        levels = np.linspace(Z.min(), Z.max(), no_of_levels)

        fig, ax = plt.subplots()

        cs = ax.contourf(X, Y, Z, levels, cmap = plt.get_cmap(cmap))
        c_bar = fig.colorbar(cs, ticks = np.append(levels[::10], levels[-1]))
        
        return fig, ax

    def get_3d_plot(self, no_of_levels = 50, cmap = 'hot'):
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        X, Y = np.meshgrid(self.grid[0], self.grid[1])
        Z = self.dispersion((X, Y))

        levels = np.linspace(Z.min(), Z.max(), no_of_levels)

        surf = ax.plot_surface(X, Y, Z, cmap=plt.get_cmap(cmap), antialiased=False)

        return fig, ax

    def get_histogram_plot(self, no_of_bins = 'auto'):
        X, Y = np.meshgrid(self.grid[0], self.grid[1])
        Z = self.dispersion((X, Y))
        Z = Z.reshape((Z.shape[0] * Z.shape[1]))

        fig, ax = plt.subplots()

        _, _, _ = ax.hist(Z, bins = no_of_bins)

        return fig, ax
