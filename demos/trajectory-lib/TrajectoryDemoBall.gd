extends KinematicBody

const LIFETIME = 5.0

export var gravity = Vector3.ZERO
export var velocity = Vector3.ZERO

var age = 0.0

func _physics_process(delta):
	age += delta
	velocity += gravity * delta
	
	if move_and_collide(velocity * delta) or age >= LIFETIME:
		queue_free()

func get_shape() -> CollisionShape:
	return $CollisionShape.shape
