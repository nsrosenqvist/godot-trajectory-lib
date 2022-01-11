#  Trajectory Lib

![Screenshot of range calculation demo](https://raw.githubusercontent.com/nsrosenqvist/godot-trajectory-lib/main/showcase/showcase-1.png)

This is a utility library that provides tools to work with, and calculate, 3D ballistic trajectories. It's a collection of helpers that are designed to cover most use-cases requiring pre-calculations of velocity, time, position or collisions.

## Credits

Much of the code is based on the trajectory solving logic of [Forrest Smith's C# library](https://github.com/forrestthewoods/lib_fts/blob/master/code/fts_ballistic_trajectory.cs) (0.1.0). Refer to [his blog](https://medium.com/@ForrestTheWoods/solving-ballistic-trajectories-b0165523348c) for extended explanations on how to use the different techniques in your project. One method is also ported over from Miziziziz's [Godot3DProjectileSolver](https://github.com/Miziziziz/Godot3DProjectileSolver/blob/master/TurretBase.gd) project.

## License

MIT

## Usage

Most techniques are explained in this [blog post](https://medium.com/@ForrestTheWoods/solving-ballistic-trajectories-b0165523348c). The code is well documented so refer to [the file](https://github.com/nsrosenqvist/godot-trajectory-lib/blob/master/addons/trajectory-lib/Trajectory.gd) for how to use the API.

Similarly to Godot's `intersect_ray`, most methods return a `Dictionary` with the relevant information of the solution. A solution object has the following keys set:

- `position`: The calculated point of impact.
- `time`: Duration of travel until the point of impact.
- `gravity`: Gravity used for solution—the methods using a fixed lateral speed return a variable gravity.
- `velocity`: The required velocity to reach the point of impact.