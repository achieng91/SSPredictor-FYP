package predictor.core.model;

import predictor.core.model.math.Vector3D;

/**
 * @author waiyan
 * BoundingPrimitives determine the shape of the Particle (Atom/Molecule/Residue)
 * Extended by BoundingBox, BoundingCube and BoundingSphere 
 */
public abstract class BoundingPrimitive {

	protected Vector3D centre;
	public abstract boolean overlap(BoundingPrimitive other);
	
	public void updateCentre(double x, double y, double z, int metric){
		centre.x = x; centre.y =y; centre.z=z; centre.metric=metric;
	}
}	