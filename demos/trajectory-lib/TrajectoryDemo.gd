extends Spatial

enum Mode {
	RANGE,
	FIXED,
	MOVING,
	FIXED_LATERAL,
	MOVING_LATERAL,
	INTERCEPT,
}

var mode = Mode.RANGE
var Ball = preload("res://demos/trajectory-lib/TrajectoryDemoBall.tscn")
var Trajectory = preload("res://addons/trajectory-lib/Trajectory.gd")

onready var ball_shape = (Ball.instance()).get_shape()
onready var _target_prev_pos = $Target.global_transform.origin
onready var _target_velocity = Vector3.ZERO

onready var player_offset = $Player.translation.y
onready var target_default = $Target.global_transform
onready var gravity_dir = ProjectSettings.get("physics/3d/default_gravity_vector")
onready var gravity_magnitude = ProjectSettings.get("physics/3d/default_gravity")
onready var gravity = gravity_dir * gravity_magnitude
onready var _default_gravity_dir = gravity_dir
onready var _default_gravity_magnitude = gravity_magnitude
onready var _default_gravity = gravity

func _unhandled_input(event):
	if event is InputEventMouseButton and event.button_index == BUTTON_LEFT and event.is_pressed():
		get_tree().set_input_as_handled()
		
		match mode:
			Mode.RANGE:
				var angle = (-$Player.global_transform.basis.z).rotated(Vector3.LEFT, deg2rad(-45.0))
				shoot(angle * $UI/Views/ProjectileSpeedHSlider.value, gravity)
			Mode.FIXED:
				var solutions = Trajectory.fixed_target($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, gravity)
				
				if solutions.empty():
					show_message()
				else:
					for s in solutions:
						shoot(s.velocity, s.gravity)
			Mode.MOVING:
				var solutions = Trajectory.moving_target($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, _target_velocity, gravity)
				
				if solutions.empty():
					show_message()
				else:
					for s in solutions:
						shoot(s.velocity, s.gravity)
			Mode.FIXED_LATERAL:
				var solution = Trajectory.fixed_target_lateral_speed($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, $UI/Views/FixedLateral/MaxHeightHSlider.value, gravity_dir)
				
				if solution.empty():
					show_message()
				else:
					shoot(solution.velocity, solution.gravity)
			Mode.MOVING_LATERAL:
				var solution = Trajectory.moving_target_lateral_speed($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, _target_velocity, $UI/Views/MovingLateral/MaxHeightHSlider.value, gravity_dir)
				
				if solution.empty():
					show_message()
				else:
					shoot(solution.velocity, solution.gravity)
			Mode.INTERCEPT:
				var pos = Trajectory.intercept_position($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, _target_velocity)
				
				if not pos:
					show_message()
				else:
					shoot($Player.global_transform.origin.direction_to(pos) * $UI/Views/ProjectileSpeedHSlider.value, Vector3.ZERO)

func show_message():
	$UI/Message.visible = true
	get_tree().create_timer(1.0).connect("timeout", $UI/Message, "set_visible", [false], CONNECT_ONESHOT)

func _ready():
	_on_Range_view_up()

func shoot(velocity : Vector3, gravity : Vector3):
	var ball = Ball.instance()
	ball.velocity = velocity
	ball.gravity = gravity
	
	add_child(ball)
	ball.global_transform.origin = $Player.global_transform.origin

func update_range():
	var initial_height = player_offset + $UI/Views/Range/InitialHeightHSlider.value
	var radius = Trajectory.range($UI/Views/ProjectileSpeedHSlider.value, gravity_magnitude, initial_height)
	
	$Player.translation.y = initial_height
	$Range.mesh.size = Vector2(radius * 2, radius * 2)
	
	var angle = (-$Player.global_transform.basis.z).rotated(Vector3.LEFT, deg2rad(-45.0))
	var velocity = angle * $UI/Views/ProjectileSpeedHSlider.value
	var samples = Trajectory.samples($Player.global_transform.origin, velocity, gravity, 3.0, 20)
	var ray_col = Trajectory.ray_collision(get_world().direct_space_state, samples, 0.08)

	if not ray_col.empty():
		var sample = samples[ray_col.sample_index]
		var ball_trans = Transform()
		ball_trans.origin = sample.position
		var col_pt = Trajectory.shape_collision(get_world().direct_space_state, ball_shape, ball_trans, sample.velocity, sample.gravity, ball_shape.radius)

		if col_pt:
			place_aim(col_pt)
			$Aim.visible = true
		else:
			place_aim(ray_col.position)
			$Aim.visible = true
	else:
		$Aim.visible = false

func place_aim(position : Vector3):
	$Aim.global_transform.origin = position + Vector3(0.0, 0.1, 0.0)

func _on_InitialHeightHSlider_value_changed(value):
	$Player.translation.y = player_offset + value
	update_range()

func _on_ProjectileSpeedHSlider_value_changed(value):
	if mode == Mode.RANGE:
		update_range()

func reset_view(selected):
	self.mode = selected

	$UI/ViewSwitcher/HBoxContainer/Range.pressed = false
	$UI/ViewSwitcher/HBoxContainer/Fixed.pressed = false
	$UI/ViewSwitcher/HBoxContainer/Moving.pressed = false
	$UI/ViewSwitcher/HBoxContainer/FixedLateral.pressed = false
	$UI/ViewSwitcher/HBoxContainer/MovingLateral.pressed = false
	$UI/ViewSwitcher/HBoxContainer/Intercept.pressed = false
	
	$UI/Views/Range.visible = false
	$UI/Views/Fixed.visible = false
	$UI/Views/Moving.visible = false
	$UI/Views/FixedLateral.visible = false
	$UI/Views/MovingLateral.visible = false
	$UI/Views/Intercept.visible = false
	
	$Aim.visible = false
	$Range.visible = false
	$UI/Intercept.visible = false
	
	$Target.visible = false
	$Target/CollisionShape.disabled = true
	$Target.global_transform = target_default
	$Target/AnimationPlayer.stop()

func _on_Range_view_up(state = true):
	reset_view(Mode.RANGE)
	$UI/ViewSwitcher/HBoxContainer/Range.pressed = true
	$UI/Views/Range.visible = true
	$Range.visible = true
	$Aim.visible = true
	update_range()

func _on_Fixed_view_up(state = true):
	reset_view(Mode.FIXED)
	$UI/ViewSwitcher/HBoxContainer/Fixed.pressed = true
	$UI/Views/Fixed.visible = true
	$Target.visible = true
	$Target/CollisionShape.disabled = false

func _on_Moving_view_up(state = true):
	reset_view(Mode.MOVING)
	$UI/ViewSwitcher/HBoxContainer/Moving.pressed = true
	$UI/Views/Moving.visible = true
	$Target.visible = true
	$Target/CollisionShape.disabled = false
	$Target/AnimationPlayer.play("Patrol")
	
func _on_FixedLateral_view_up(state = true):
	reset_view(Mode.FIXED_LATERAL)
	$UI/ViewSwitcher/HBoxContainer/FixedLateral.pressed = true
	$UI/Views/FixedLateral.visible = true
	$Target.visible = true
	$Target/CollisionShape.disabled = false

func _on_MovingLateral_view_up(state = true):
	reset_view(Mode.MOVING_LATERAL)
	$UI/ViewSwitcher/HBoxContainer/MovingLateral.pressed = true
	$UI/Views/MovingLateral.visible = true
	$Target.visible = true
	$Target/CollisionShape.disabled = false
	$Target/AnimationPlayer.play("Patrol")

func _on_Intercept_view_up(state = true):
	reset_view(Mode.INTERCEPT)
	$UI/ViewSwitcher/HBoxContainer/Intercept.pressed = true
	$UI/Views/Intercept.visible = true
	$Target.visible = true
	$Target/CollisionShape.disabled = false
	$Target/AnimationPlayer.play("Patrol")
	$UI/Intercept.visible = true
	
func _physics_process(delta):
	_target_velocity = ($Target.global_transform.origin - _target_prev_pos) / delta
	_target_prev_pos = $Target.global_transform.origin
	
	if mode == Mode.INTERCEPT and $UI/Intercept/VisibilityNotifier.is_on_screen() and _target_velocity.length() > 0.0:
		var intercept_time = Trajectory.intercept_time($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, _target_velocity)
		var intercept_pos = Trajectory.intercept_position($Player.global_transform.origin, $UI/Views/ProjectileSpeedHSlider.value, $Target.global_transform.origin, _target_velocity)
		
		if intercept_pos:
			$UI/Intercept.visible = true
			$UI/Intercept.position = $Camera.unproject_position(intercept_pos)
			$UI/Intercept/TimeLabel.text = str(intercept_time)
		else:
			$UI/Intercept.visible = false
	else:
		$UI/Intercept.visible = false

func _on_GravitySwitcher_toggled(button_pressed):
	if button_pressed:
		if mode == Mode.RANGE:
			_on_Fixed_view_up(true)
		
		$GravityField/CollisionShape.disabled = false
		$UI/ViewSwitcher/HBoxContainer/Range.disabled = true
		gravity_dir = $GravityField.gravity_vec
		gravity_magnitude = $GravityField.gravity
		gravity = gravity_dir * gravity_magnitude
	else:
		$GravityField/CollisionShape.disabled = true
		$UI/ViewSwitcher/HBoxContainer/Range.disabled = false
		gravity = _default_gravity
		gravity_dir = _default_gravity_dir
		gravity_magnitude = _default_gravity_magnitude
