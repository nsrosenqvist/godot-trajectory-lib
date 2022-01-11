extends Node

# Much of the ballistic trajectory solving code is ported over from Forrest Smith's
# public domain C# library (version 0.1.0):
# 	https://github.com/forrestthewoods/lib_fts/blob/master/code/fts_ballistic_trajectory.cs
#
# These methods includes:
# 	- get_cubic_root [GetCubicRoot]
# 	- solve_quadric [SolveQuadric]
# 	- solve_cubic [SolveCubic]
# 	- solve_quartic [SolveQuartic]
# 	- range [ballistic_range]
# 	- fixed_target [solve_ballistic_arc:258]
# 	- moving_target [solve_ballistic_arc:354]
# 	- fixed_target_lateral_speed [solve_ballistic_arc_lateral:454]
# 	- moving_target_lateral_speed [solve_ballistic_arc_lateral:501]
#
# Additional reading that explains each function and its use cases can be found at:
# 	https://www.forrestthewoods.com/blog/solving_ballistic_trajectories/
#
# The methods solve_quadric, solve_cubic and solve_quartic were originally written
# in C by Jochen Schwarze (schwarze@isa.de) and ported over to C# by Forrest Smith:
# 	https://github.com/erich666/GraphicsGems/blob/240a34f2ad3fa577ef57be74920db6c4b00605e4/gems/Roots3And4.c

# intercept_time is adapted from Miziziziz's calculate_direction of
# Godot3DProjectileSolver (MIT Copyright (c) 2021 Miziziziz). This implementation
# also fixes a bug present in original code (see line 47 and 49):
# 	https://github.com/Miziziziz/Godot3DProjectileSolver/blob/master/TurretBase.gd

## 
## A library providing functions to help calculating and working with 3D
## projectile trajectories
## 
## @desc:
##     This is a collection of helper methods that is designed to cover most
##     use-cases requiring pre-calculating, or testing, 3D projectile trajectories.
##     The code is MIT-licensed and is based primarily on the C# ballistic
##     trajectory library created by Forrest Smith.
## 
## @tutorial: https://www.github.com/nsrosenqvist/godot-trajectory-lib#readme
##

## Retrieve cubic root
static func get_cubic_root(value : float) -> float:
	if value > 0.0:
		return pow(value, 1.0 / 3.0)
	elif value < 0:
		return -pow(-value, 1.0 / 3.0)
	else:
		return 0.0

## Solve quadratic equation: c0*x^2 + c1*x + c2. 
## Returns solutions as an [Array].
static func solve_quadric(c0 : float, c1 : float, c2 : float) -> Array:
	var p : float
	var q : float
	var D : float
	
	# normal form: x^2 + px + q = 0
	p = c1 / (2 * c0);
	q = c2 / c0;

	D = p * p - q;

	if is_zero_approx(D):
		return [-p]
	elif D < 0:
		return []
	else: # elif D > 0:
		var sqrt_D    : float = sqrt(D)
		var solutions : Array = []
		
		solutions.append(sqrt_D - p)
		solutions.append(-sqrt_D - p)
		return solutions

## Solve cubic equation: c0*x^3 + c1*x^2 + c2*x + c3. 
## Returns solutions as an [Array].
static func solve_cubic(c0 : float, c1 : float, c2 : float, c3 : float) -> Array:
	var solutions : Array = []

	var sub : float
	var A : float
	var B : float
	var C : float
	var sq_A : float
	var p : float
	var q : float
	var cb_p : float
	var D : float

	# normal form: x^3 + Ax^2 + Bx + C = 0
	A = c1 / c0
	B = c2 / c0
	C = c3 / c0

	# substitute x = y - A/3 to eliminate quadric term:  x^3 +px + q = 0
	sq_A = A * A
	p = 1.0/3 * (- 1.0/3 * sq_A + B)
	q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C)

	# use Cardano's formula
	cb_p = p * p * p
	D = q * q + cb_p

	if is_zero_approx(D):
		if is_zero_approx(q): # one triple solution
			solutions.append(0)
		else: # one single and one double solution
			var u : float = get_cubic_root(-q)
			solutions.append(2 * u)
			solutions.append(-u)
	elif D < 0: # Casus irreducibilis: three real solutions
		var phi : float = 1.0/3 * acos(-q / sqrt(-cb_p))
		var t   : float = 2 * sqrt(-p)

		solutions.append(t * cos(phi))
		solutions.append(-t * cos(phi + PI / 3))
		solutions.append(-t * cos(phi - PI / 3))
	else: # one real solution
		var sqrt_D : float = sqrt(D)
		var u      : float = get_cubic_root(sqrt_D - q)
		var v      : float = -get_cubic_root(sqrt_D + q)

		solutions.append(u + v)

	# resubstitute
	sub = 1.0/3 * A

	for i in range(0, solutions.size()):
		solutions[i] -= sub

	return solutions

## Solve quartic function: c0*x^4 + c1*x^3 + c2*x^2 + c3*x + c4. 
## Returns solutions as an [Array].
static func solve_quartic(c0 : float, c1 : float, c2 : float, c3 : float, c4 : float) -> Array: #out double s0, out double s1, out double s2, out double s3)
	var solutions : Array = []

	# TODO: Use typed array when GDScript 2.0 releases
	# https://github.com/godotengine/godot-proposals/issues/2181
	# var coeffs : Array[float] = []
	var coeffs : Array = [NAN, NAN, NAN, NAN]
	var z : float
	var u : float
	var v : float
	var sub : float
	var A : float
	var B : float
	var C : float
	var D : float
	var sq_A : float
	var p : float
	var q : float
	var r : float

	# normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0
	A = c1 / c0
	B = c2 / c0
	C = c3 / c0
	D = c4 / c0

	# substitute x = y - A/4 to eliminate cubic term: x^4 + px^2 + qx + r = 0
	sq_A = A * A
	p = - 3.0/8 * sq_A + B
	q = 1.0/8 * sq_A * A - 1.0/2 * A * B + C
	r = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D

	if is_zero_approx(r):
		# no absolute term: y(y^3 + py + q) = 0
		coeffs[3] = q
		coeffs[2] = p
		coeffs[1] = 0
		coeffs[0] = 1
		
		solutions = solve_cubic(coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	else:
		# solve the resolvent cubic ...
		coeffs[3] = 1.0/2 * r * p - 1.0/8 * q * q
		coeffs[2] = - r
		coeffs[1] = - 1.0/2 * p
		coeffs[0] = 1

		solutions = solve_cubic(coeffs[0], coeffs[1], coeffs[2], coeffs[3])

		# ... and take the one real solution ...
		z = solutions.front()

		# ... to build two quadric equations
		u = z * z - r
		v = 2 * z - p

		if is_zero_approx(u):
			u = 0
		elif u > 0:
			u = sqrt(u)
		else:
			return []

		if is_zero_approx(v):
			v = 0
		elif v > 0:
			v = sqrt(v)
		else:
			return []

		coeffs[2] = z - u
		coeffs[1] = -v if q < 0 else v
		coeffs[0] = 1

		solutions = solve_quadric(coeffs[0], coeffs[1], coeffs[2])

		coeffs[2]= z + u
		coeffs[1] = v if q < 0 else -v
		coeffs[0] = 1
		
		solutions.append_array(solve_quadric(coeffs[0], coeffs[1], coeffs[2]))

	# resubstitute
	sub = 1.0/4 * A

	for i in range(0, solutions.size()):
		solutions[i] -= sub

	return solutions

## Calculate the maximum range that a ballistic projectile can be fired on given speed and gravity.
##
## speed [float]: projectile velocity
## gravity_magnitude [float]: force of gravity
## initial_height [float]: distance above flat terrain
##
## return [float]: maximum range
static func range(speed : float, gravity_magnitude : float, initial_height : float) -> float:
	assert(speed > 0 and gravity_magnitude > 0.0, "range called with invalid data.");

	# Derivation
	#	(1) x = speed * time * cos O
	#	(2) y = initial_height + (speed * time * sin O) - (.5 * gravity*time*time)
	#	(3) via quadratic: t = (speed*sin O)/gravity + sqrt(speed*speed*sin O + 2*gravity*initial_height)/gravity    [ignore smaller root]
	#	(4) solution: range = x = (speed*cos O)/gravity * sqrt(speed*speed*sin O + 2*gravity*initial_height)    [plug t back into x=speed*time*cos O]
	var angle : float = deg2rad(45.0) # no air resistence, so 45 degrees provides maximum range
	var a_cos : float = cos(angle)
	var a_sin : float = sin(angle)

	return (speed*a_cos/gravity_magnitude) * (speed*a_sin + sqrt(speed*speed*a_sin*a_sin + 2*gravity_magnitude*initial_height))

## Solve firing angles for a ballistic projectile with speed and gravity to hit a fixed position.
##
## proj_pos [Vector3]: point projectile will fire from
## proj_speed [float]: scalar speed of projectile
## target_pos [Vector3]: point projectile is trying to hit
## gravity [Vector3]: gravity vector
##
## Returns 0-2 solutions ([Array]) in ascending order of angle incline, low to
## high. Each solution is a [Dictionary] with velocity, gravity, position, and time.
static func fixed_target(proj_pos : Vector3, proj_speed : float, target_pos : Vector3, gravity : Vector3) -> Array:
	# Support calculations relative to gravity
	var gravity_dir       : Vector3 = gravity.normalized()
	var gravity_magnitude : float   = gravity.length()
	var gravity_basis     : Basis

	if gravity_dir != Vector3.DOWN:
		gravity_basis = _gravity_basis(gravity_dir)
		proj_pos = gravity_basis.xform_inv(proj_pos)
		target_pos = gravity_basis.xform_inv(target_pos)

	assert(proj_pos != target_pos and proj_speed > 0.0 and gravity_magnitude > 0.0, "fixed_target called with invalid data.");

	# Derivation
	#	(1) x = v*t*cos O
	#	(2) y = v*t*sin O - .5*g*t^2
	#
	#	(3) t = x/(cos O*v)                                        [solve t from (1)]
	#	(4) y = v*x*sin O/(cos O * v) - .5*g*x^2/(cos^2 O*v^2)     [plug t into y=...]
	#	(5) y = x*tan O - g*x^2/(2*v^2*cos^2 O)                    [reduce; cos/sin = tan]
	#	(6) y = x*tan O - (g*x^2/(2*v^2))*(1+tan^2 O)              [reduce; 1+tan O = 1/cos^2 O]
	#	(7) 0 = ((-g*x^2)/(2*v^2))*tan^2 O + x*tan O - (g*x^2)/(2*v^2) - y    [re-arrange]
	#	Quadratic! a*p^2 + b*p + c where p = tan O
	#
	#	(8) let gxv = -g*x*x/(2*v*v)
	#	(9) p = (-x +- sqrt(x*x - 4gxv*(gxv - y)))/2*gxv           [quadratic formula]
	#	(10) p = (v^2 +- sqrt(v^4 - g(g*x^2 + 2*y*v^2)))/gx        [multiply top/bottom by -2*v*v/x; move 4*v^4/x^2 into root]
	#	(11) O = atan(p)
	
	var diff    : Vector3 = target_pos - proj_pos
	var diff_xz : Vector3 = Vector3(diff.x, 0.0, diff.z)
	var ground_dist = diff_xz.length()

	var speed2 : float = proj_speed*proj_speed
	var speed4 : float = proj_speed*proj_speed*proj_speed*proj_speed
	var y      : float = diff.y
	var x      : float = ground_dist
	var gx     : float = gravity_magnitude*x

	var root : float = speed4 - gravity_magnitude*(gravity_magnitude*x*x + 2*y*speed2)

	# No solution
	if root < 0:
		return []

	root = sqrt(root)

	var solutions : Array = []
	var low_ang    : float = atan2(speed2 - root, gx)
	var high_ang   : float = atan2(speed2 + root, gx)

	var ground_dir : Vector3 = diff_xz.normalized()
	var low_velocity = ground_dir*cos(low_ang)*proj_speed + Vector3.UP*sin(low_ang)*proj_speed
	
	solutions.append({
		"position": target_pos,
		"time": ground_dist / (cos(low_ang) * proj_speed),
		"velocity": low_velocity,
		"gravity": gravity,
	})

	if low_ang != high_ang:
		var high_velocity = ground_dir*cos(high_ang)*proj_speed + Vector3.UP*sin(high_ang)*proj_speed
		
		solutions.append({
			"position": target_pos,
			"time": ground_dist / (cos(high_ang) * proj_speed),
			"velocity": high_velocity,
			"gravity": gravity,
		})

	# Make solutions global again
	if gravity_dir != Vector3.DOWN:
		for i in range(0, solutions.size()):
			solutions[i].position = gravity_basis.xform(solutions[i].position)
			solutions[i].velocity = gravity_basis.xform(solutions[i].velocity)
	
	return solutions

static func _gravity_basis(gravity_dir : Vector3) -> Basis:
	var basis    := Basis()
	var to_cross := Vector3.LEFT
	
	if to_cross == gravity_dir:
		to_cross = Vector3.FORWARD
	
	basis.y = -(gravity_dir)
	basis.x = basis.y.cross(to_cross)
	basis.z = basis.y.cross(basis.x)
	
	return basis.orthonormalized()

## Solve firing angles for a ballistic projectile with speed and gravity to hit a target moving with
## constant, linear velocity.
##
## proj_pos [Vector3]: point projectile will fire from
## proj_speed [float]: scalar speed of projectile
## target_pos [Vector3]: point projectile is trying to hit
## target_pos [Vector3]: current velocity of target
## gravity [Vector3]: gravity vector
##
## Returns 0-2 solutions ([Array]) in ascending order of angle incline, low to
## high. Each solution is a [Dictionary] with velocity, gravity, position, and time.
static func moving_target(proj_pos : Vector3, proj_speed : float, target_pos : Vector3, target_velocity : Vector3, gravity : Vector3) -> Array:
	# Support calculations relative to gravity
	var gravity_dir       : Vector3 = gravity.normalized()
	var gravity_magnitude : float   = gravity.length()
	var gravity_basis     : Basis
	
	if gravity_dir != Vector3.DOWN:
		gravity_basis = _gravity_basis(gravity_dir)
		proj_pos = gravity_basis.xform_inv(proj_pos)
		target_pos = gravity_basis.xform_inv(target_pos)
		target_velocity = gravity_basis.xform_inv(target_velocity)

	# Derivation 
	#
	# For full derivation see: blog.forrestthewoods.com
	# Here is an abbreviated version.
	#
	# Four equations, four unknowns (solution.x, solution.y, solution.z, time):
	#	(1) proj_pos.x + solution.x*time = target_pos.x + target_vel.x*time
	#	(2) proj_pos.y + solution.y*time + .5*G*t = target_pos.y + target_vel.y*time
	#	(3) proj_pos.z + solution.z*time = target_pos.z + target_vel.z*time
	#	(4) proj_speed^2 = solution.x^2 + solution.y^2 + solution.z^2
	#
	#	(5) Solve for solution.x and solution.z in equations (1) and (3)
	#	(6) Square solution.x and solution.z from (5)
	#	(7) Solve solution.y^2 by plugging (6) into (4)
	#	(8) Solve solution.y by rearranging (2)
	#	(9) Square (8)
	#	(10) Set (8) = (7). All solution.xyz terms should be gone. Only time remains.
	#	(11) Rearrange 10. It will be of the form a*^4 + b*t^3 + c*t^2 + d*t * e. This is a quartic.
	#	(12) Solve the quartic using SolveQuartic.
	#	(13) If there are no positive, real roots there is no solution.
	#	(14) Each positive, real root is one valid solution
	#	(15) Plug each time value into (1) (2) and (3) to calculate solution.xyz
	#	(16) The end.

	var G : float = gravity_magnitude

	var A : float = proj_pos.x
	var B : float = proj_pos.y
	var C : float = proj_pos.z
	var M : float = target_pos.x
	var N : float = target_pos.y
	var O : float = target_pos.z
	var P : float = target_velocity.x
	var Q : float = target_velocity.y
	var R : float = target_velocity.z
	var S : float = proj_speed

	var H : float = M - A
	var J : float = O - C
	var K : float = N - B
	var L : float = -0.5 * G

	# Quartic Coeffecients
	var c0 : float = L*L
	var c1 : float = -2*Q*L
	var c2 : float = Q*Q - 2*K*L - S*S + P*P + R*R
	var c3 : float = 2*K*Q + 2*H*P + 2*J*R
	var c4 : float = K*K + H*H + J*J

	# Solve quartic
	# TODO: Use typed array when GDScript 2.0 releases
	# https://github.com/godotengine/godot-proposals/issues/2181
	# var times : Array[float] = []
	var times : Array = solve_quartic(c0, c1, c2, c3, c4)
	
	# Sort so faster collision is found first
	times.sort()

	# Plug quartic solutions into base equations
	# There should never be more than 2 positive, real roots.
	# TODO: Use typed array when GDScript 2.0 releases
	# https://github.com/godotengine/godot-proposals/issues/2181
	# var solutions : Array[Vector3] = []
	var solutions : Array = []

	for t in times:
		if t <= 0 or is_nan(t):
			continue

		var velocity : Vector3 = Vector3.ZERO
		velocity.x = float((H+P*t)/t)
		velocity.y = float((K+Q*t-L*t*t)/ t)
		velocity.z = float((J+R*t)/t)
		
		solutions.append({
			"position": position_at(t, proj_pos, velocity, gravity),
			"time": t,
			"velocity": velocity,
			"gravity": gravity,
		})
		
		if solutions.size() >= 2:
			break

	# Make out parameters global again
	if gravity_dir != Vector3.DOWN:
		for i in range(0, solutions.size()):
			solutions[i].velocity = gravity_basis.xform(solutions[i].velocity)
			solutions[i].position = gravity_basis.xform(solutions[i].position)

	return solutions

## Solve the firing arc with a fixed lateral speed. Vertical speed and gravity varies. 
## This enables a visually pleasing arc.
##
## proj_pos [Vector3]: point projectile will fire from
## lateral_speed [float]: scalar speed of projectile along XZ plane
## target_pos [Vector3]: point projectile is trying to hit
## max_height [float]: height above max(proj_pos, impact_pos) for projectile to peak at
## gravity_dir [Vector3]: direction of gravity
##
## Returns a [Dictionary] with velocity, gravity, position, and time, or empty if
## no solution was found.
static func fixed_target_lateral_speed(proj_pos : Vector3, lateral_speed : float, target_pos : Vector3, max_height : float, gravity_dir : Vector3) -> Dictionary:
	# Support calculations relative to gravity
	var gravity_basis : Basis
	
	if gravity_dir != Vector3.DOWN:
		gravity_basis = _gravity_basis(gravity_dir)
		proj_pos = gravity_basis.xform_inv(proj_pos)
		target_pos = gravity_basis.xform_inv(target_pos)
	
	assert(proj_pos != target_pos and lateral_speed > 0 and max_height > proj_pos.y, "fixed_target_lateral_speed called with invalid data.");

	var solution : Dictionary = {
		"velocity": Vector3.ZERO,
		"gravity": Vector3.ZERO,
		"position": target_pos,
		"time": 0.0,
	}
	
	var diff         : Vector3 = target_pos - proj_pos
	var diff_xz      : Vector3 = Vector3(diff.x, 0.0, diff.z)
	var lateral_dist : float   = diff_xz.length()

	if lateral_dist == 0.0:
		return {}

	var time : float = lateral_dist / lateral_speed;

	solution.velocity = diff_xz.normalized() * lateral_speed

	# System of equations. Hit max_height at t=.5*time. Hit target at t=time.
	#
	# peak = y0 + vertical_speed*halfTime + .5*gravity*halfTime^2
	# end = y0 + vertical_speed*time + .5*gravity*time^s
	# Wolfram Alpha: solve b = a + .5*v*t + .5*g*(.5*t)^2, c = a + vt + .5*g*t^2 for g, v
	var a : float = proj_pos.y # initial
	var b : float = max_height # peak
	var c : float = target_pos.y # final

	var gravity_magnitude = -4*(a - 2*b + c) / (time*time)
	
	solution.velocity.y = -(3*a - 4*b + c) / time
	solution.time = time
	solution.gravity = Vector3.DOWN * gravity_magnitude
	
	# Make solution global again
	if gravity_dir != Vector3.DOWN:
		solution.gravity = gravity_basis.xform(solution.gravity)
		solution.velocity = gravity_basis.xform(solution.velocity)

	return solution

## Solve the firing arc with a fixed lateral speed. Vertical speed and gravity varies. 
## This enables a visually pleasing arc.
##
## proj_pos [Vector3]: point projectile will fire from
## lateral_speed [float]: scalar speed of projectile along XZ plane
## target_pos [Vector3]: point projectile is trying to hit
## max_height [float]: height above Max(proj_pos, impact_pos) for projectile to peak at
##
## Returns a [Dictionary] with velocity, gravity, position, and time, or empty if
## no solution was found.
static func moving_target_lateral_speed(proj_pos : Vector3, lateral_speed : float, target_pos : Vector3, target_velocity : Vector3, max_height : float, gravity_dir : Vector3) -> Dictionary:
	# Support calculations relative to gravity
	var gravity_basis : Basis
	
	if gravity_dir != Vector3.DOWN:
		gravity_basis = _gravity_basis(gravity_dir)
		proj_pos = gravity_basis.xform_inv(proj_pos)
		target_pos = gravity_basis.xform_inv(target_pos)
		target_velocity = gravity_basis.xform_inv(target_velocity)
	
	assert(proj_pos != target_pos and lateral_speed > 0.0, "moving_target_lateral_speed called with invalid data.")

	# Initialize output variables
	var solution : Dictionary = {
		"velocity": Vector3.ZERO,
		"gravity": 0.0,
		"position": Vector3.ZERO,
		"time": 0.0,
	}

	# Ground plane terms
	var target_vel_xz : Vector3 = Vector3(target_velocity.x, 0.0, target_velocity.z)
	var diff_xz       : Vector3 = target_pos - proj_pos
	diff_xz.y = 0.0

	# Derivation
	#	(1) Base formula: |P + V*t| = S*t
	#	(2) Substitute variables: |diffXZ + targetVelXZ*t| = S*t
	#	(3) Square both sides: Dot(diffXZ,diffXZ) + 2*Dot(diffXZ, targetVelXZ)*t + Dot(targetVelXZ, targetVelXZ)*t^2 = S^2 * t^2
	#	(4) Quadratic: (Dot(targetVelXZ,targetVelXZ) - S^2)t^2 + (2*Dot(diffXZ, targetVelXZ))*t + Dot(diffXZ, diffXZ) = 0
	var c0 : float = target_vel_xz.dot(target_vel_xz) - lateral_speed*lateral_speed
	var c1 : float = 2.0 * diff_xz.dot(target_vel_xz)
	var c2 : float = diff_xz.dot(diff_xz)
	
	var solutions : Array = solve_quadric(c0, c1, c2)
	var n  : int   = solutions.size()
	var t0 : float = solutions[0] if n > 0 else NAN
	var t1 : float = solutions[1] if n > 1 else NAN

	# pick smallest, positive time
	var valid0 : bool = (n > 0 and t0 > 0)
	var valid1 : bool = (n > 1 and t1 > 0)

	var t : float
	
	if not valid0 and not valid1:
		return {}
	elif valid0 and valid1:
		t = min(float(t0), float(t1))
	else:
		t = float(t0) if valid0 else float(t1)

	# Calculate impact point
	solution.position = target_pos + (target_velocity*t)

	# Calculate fire velocity along XZ plane
	var dir : Vector3 = solution.position - proj_pos
	solution.velocity = Vector3(dir.x, 0.0, dir.z).normalized() * lateral_speed

	# Solve system of equations. Hit max_height at t=.5*time. Hit target at t=time.
	#
	# peak = y0 + vertical_speed*halfTime + .5*gravity*halfTime^2
	# end = y0 + vertical_speed*time + .5*gravity*time^s
	# Wolfram Alpha: solve b = a + .5*v*t + .5*g*(.5*t)^2, c = a + vt + .5*g*t^2 for g, v
	var a : float = proj_pos.y # initial
	var b : float = max(proj_pos.y, solution.position.y) + max_height # peak
	var c : float = solution.position.y # final

	var gravity_magnitude : float = -4*(a - 2*b + c) / (t*t)
	
	solution.velocity.y = -(3*a - 4*b + c) / t
	solution.gravity = Vector3.DOWN * gravity_magnitude
	solution.time = t
	
	# Make solution global again
	if gravity_dir != Vector3.DOWN:
		solution.gravity = gravity_basis.xform(solution.gravity)
		solution.position = gravity_basis.xform(solution.position)
		solution.velocity = gravity_basis.xform(solution.velocity)

	return solution

## Calculates the time required to intercept another object
## 
## origin [Vector3]: Position of follower
## velocity [Vector3]: Velocity of follower
## target [Vector3]: Position of target
## target_velocity [Vector3]: Velocity of target
##
## Returns either the time ([float]) or NAN.
static func intercept_time(origin : Vector3, speed : float, target : Vector3, target_velocity := Vector3.ZERO) -> float:
	if target_velocity == Vector3.ZERO:
		return origin.distance_to(target) / speed if speed > 0.0 else NAN
	else:
		var Pti : Vector3 = target
		var Pbi : Vector3 = origin
		var D   : float   = Pti.distance_to(Pbi)
		var Vt  : Vector3 = target_velocity
		var St  : float   = Vt.length()
		var Sb  : float   = speed
		
		var cos_theta : float = Pti.direction_to(Pbi).dot(Vt.normalized())
		var root      : float = sqrt(2*D*St*cos_theta + 4*(Sb*Sb - St*St)*D*D )
		var t1        : float = (-2*D*St*cos_theta + root) / (2*(Sb*Sb - St*St))
		var t2        : float = (-2*D*St*cos_theta - root) / (2*(Sb*Sb - St*St))
		var t         : float = min(t1, t2)
		
		if t < 0:
			t = max(t1, t2)
		if t < 0:
			return NAN # can't intercept, target too fast

		return t

## Calculates the position of an object at the time it intercepts its target
## 
## origin [Vector3]: Position of follower
## velocity [Vector3]: Velocity of follower
## target [Vector3]: Position of target
## target_velocity [Vector3]: Velocity of target
##
## Returns either the position ([Vector3]) or null.
static func intercept_position(origin : Vector3, speed : float, target : Vector3, target_velocity := Vector3.ZERO):
	if target_velocity == Vector3.ZERO:
		return target
	
	var time : float = intercept_time(origin, speed, target, target_velocity)
	
	if not is_nan(time):
		return target + (target_velocity * time)

	return null

static func _float_range(from : float, to = null, step = 1.0) -> Array:
	if to == null:
		to = from
		from = 0.0
	
	var result  : Array = []
	var current : float = from
	
	while current <= (to - step):
		result.append(current)
		current += step
	
	return result

## Calculates sample points in space along an object's trajectory.
## 
## origin [Vector3]: Start of trajectory
## velocity [Vector3]: Velocity of object
## gravity [Vector3]: Gravity acceleration
## time [float]: Duration of travel that the curve will cover
## iterations [int]: Number of sample points that will be calculated
##
## Returns an [Array] of [Dictionary]'s that specify:
##     - position (in space)
##     - time (along trajectory)
##     - velocity (at time)
##     - gravity (constant for calculation)
static func samples(origin : Vector3, velocity : Vector3, gravity : Vector3, time : float, iterations : int = 10) -> Array:
	var samples   : Array = []
	var time_step : float = time / iterations

	var s0 = origin
	var v0 = velocity
	var a = gravity

	for t in _float_range(0.0, time + time_step, time_step):
		var s = s0 + v0*t + 0.5*a*pow(t, 2.0)
		var v = v0 + a*t
		
		samples.append({
			"position": s,
			"time": t,
			"velocity": v,
			"gravity": a,
		})

	return samples

## Get an estimated position after the specified set amount of time
## providing that no collisions have interrupted the trajectory.
## 
## origin [Vector3]: Start of trajectory
## time [float]: Time of position to solve for
## velocity [Vector3]: Velocity of object
## gravity [Vector3]: Gravity acceleration
##
## Returns the position as a [Vector3]
static func position_at(time : float, origin : Vector3, velocity : Vector3, gravity : Vector3) -> Vector3:
	return origin + velocity*time + 0.5*gravity*pow(time, 2.0)

## Get an estimated velocity after the specified set amount of time
## providing that no collisions have interrupted the trajectory.
## 
## time [float]: Time of velocity to solve for
## velocity [Vector3]: Velocity of object
## gravity [Vector3]: Gravity acceleration
##
## Returns the velocity as a [Vector3]
static func velocity_at(time : float, velocity : Vector3, gravity : Vector3) -> Vector3:
	return velocity + gravity * time

## Tests when a ray would collide along the provided trajectory. The function
## casts a ray in-between each set of sample points along the trajectory and 
## should therefore be considered a rough estimation, with its accuracy 
## depending on how many samples the trajectory was created with.
## 
## space_state [PhysicsDirectSpaceState]: The physics space that will be tested
## samples [Array]: Trajectory point samples to test (see trajectory_samples)
## collision_margin [float]: Helps make sure that the ray actually hits if a
##     sample point is precisely on the edge of a physical object.
## exclude [Array]: List of nodes to exclude when casting rays
## collision_mask [int]: The bitmap of the ray's collision_mask
## collide_with_bodies [bool]: Flag to indicate whether to test against physics bodies
## collide_with_areas [bool]: Flag to indicate whether to test against physics areas
##
## Returns the collision as a [Dictionary]; empty if the ray never hit anything.
## In addition to the normal collision info, the index of the colliding raycast
## is provided as well as "sample_index".
static func ray_collision(space_state : PhysicsDirectSpaceState, samples : Array, collision_margin : float = 0.04, exclude : Array = [], collision_mask: int = 0x7FFFFFFF, collide_with_bodies : bool = true, collide_with_areas : bool = false) -> Dictionary:
	var num : int        = samples.size()
	var col : Dictionary = {}
	
	assert(samples.size() > 1, "ray_collision: Sample size must be higher than one.")
	
	for i in range(0, num - 1):
		var p = samples[i].position
		var n = samples[i + 1].position
		
		var dir = (n - p).normalized()
		p += (-dir) * collision_margin
		n += dir * collision_margin
		
		col = space_state.intersect_ray(p, n, exclude, collision_mask, collide_with_bodies, collide_with_areas)
		
		if not col.empty():
			col["sample_index"] = i
			break
	
	return col

## Tests when a shape would collide along its trajectory. Since this method is
## more resource-heavy than checking for collisions using a raycast its recommended
## to first retrieve the colliding index using trajectory_ray_collision and only
## perform the simulation from that point in space.
## 
## space_state [PhysicsDirectSpaceState]: The physics space that will be tested
## shape [Shape]: The physics shape that will be used for collisions
## origin [Vector3]: The start position of the trajectory
## velocity [Vector3]: Starting velocity of object
## gravity [Vector3]: Gravity affecting the object
## collider_radius [float]: Collision point offset. If zero, then the origin of
##     the collider will be returned and not the actual point of impact.
## time_step [float]: Time interval that will be tested, by default the value is
##     the project defined physics fps
## max_time [float]: Time limit for simulation.
## exclude [Array]: List of nodes to exclude when casting rays
## collision_mask [int]: The bitmap of the ray's collision_mask
## collide_with_bodies [bool]: Flag to indicate whether to test against physics bodies
## collide_with_areas [bool]: Flag to indicate whether to test against physics areas
##
## Returns the collision point as a [Vector3], or null if no collision.
static func shape_collision(space_state : PhysicsDirectSpaceState, shape : Shape, transform : Transform, velocity : Vector3, gravity : Vector3, collider_radius : float = 0.0, max_time : float = 10.0, time_step : float = 1.0 / ProjectSettings.get_setting("physics/common/physics_fps"), exclude : Array = [], collision_mask: int = 0x7FFFFFFF, collide_with_bodies : bool = true, collide_with_areas : bool = false, num_collisions : int = 1):
	assert(max_time > 0.0, "trajectory_shape_collision max_time must be a positive value.")

	var params : PhysicsShapeQueryParameters = PhysicsShapeQueryParameters.new()
	params.exclude = exclude
	params.collision_mask = collision_mask
	params.collide_with_bodies = collide_with_bodies
	params.collide_with_areas = collide_with_areas
	params.set_shape(shape)
	params.set_transform(transform)
	
	var current = 0.0
	var origin = transform.origin
	
	while current <= max_time:
		transform.origin = origin
		params.set_transform(transform)
		var motion = space_state.cast_motion(params, velocity * time_step)

		if motion[1] < 1.0:
			return transform.origin + (velocity * time_step * motion[1]) + (velocity.normalized() * collider_radius)
		
		velocity += gravity * time_step
		origin += velocity * time_step
	
	return null
