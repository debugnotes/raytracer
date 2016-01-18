using System;
using System.Collections;
using System.Drawing;
using System.Text;
using MatrixUtils;

namespace Raytracer
{
	class Raytracer
	{
		private ArrayList arrObjects = null;
		private ArrayList arrPointLightSources = null;
		private Camera camEye = null;
		private int iMaxRecursion = 1;
		private Color clBackground = Color.Black;
		private Vector vectAmbientLightIntensity = new Vector(0.5, 0.5, 0.5); 
		private Vector vectAmbientLightColor = new Vector(1.0, 1.0, 1.0);  

		public Raytracer()
		{
			camEye = new Camera(0, 0, 200, 0, 0, 0, 0, 1, 0);
		}

		public Camera Eye
		{
			get
			{
				return camEye;
			}
		}
		/// <summary>
		/// Renders a scene using ray tracing.
		/// </summary>

		///
		public Bitmap Render(int iWidth, int iHeight)
		{
			// setup a bitmap object to receive the drawing commands
			Bitmap bmpBitmap = new Bitmap(iWidth, iHeight,
				System.Drawing.Imaging.PixelFormat.Format24bppRgb);

			// setup the camera to have a 45 degree view angle
			// in the X and Y plane
			double dTheta = 45.0 * Math.PI / 180.0 / 2.0; // 45/2 degrees, in radians
			double dLeft = Math.Tan(dTheta) * (-iWidth/2);
			double dRight = Math.Tan(dTheta) * (iWidth / 2);
			double dTop = Math.Tan(dTheta) * (iHeight / 2);
			double dBottom = Math.Tan(dTheta) * (-iHeight / 2);
			double dNear = 0;
			//double dFar = -1000;
			
			// get the U, V and W basis vectors from the Camera
			Vector vectU = camEye.U;
			Vector vectV = camEye.V;
			Vector vectW = camEye.W;

			for (int iY = 0; iY < iHeight; iY++)
			{
				double dvs = dTop + ((dBottom - dTop) * ((iY + 0.5) / iHeight));
				for (int iX = 0; iX < iWidth; iX++)
				{
					double dus = dLeft + ((dRight - dLeft) * ((iX + 0.5) / iWidth));
					double dws = dNear;
					// fire a primary ray from the camera onto the near plane of the 
					// viewing window, adjusting the UVW coordinates of the ray to world
					// coordinates
					Vector vectRayEnd = ((dus * vectU) + (dvs * vectV) + (dws * vectW)); // -camEye.Position;
						//;
					Vector vectRayStart = camEye.Position;
					// convert the ray from camera coordinates to world coordinates
					bmpBitmap.SetPixel(iX, iY, VectorToColor(TraceRay(vectRayStart,
						vectRayEnd, 0)));
				}
			}
			return bmpBitmap;
		}

		private Vector TraceRay(Vector vectRayStart, Vector vectRayEnd, int iRayCount)
		{
			// establish a vector of base intenisites (current "black").
			Vector vectIntensity = new Vector(3);
			Vector vectSpecular = new Vector(3);
			// make sure we don't recurse any more than we have to
			if (iRayCount > iMaxRecursion)
				return vectIntensity;
			// first, find the object that intersects the ray at the start position
			Object3D obj3dClosest = null;
			double dtClosest = 0;
			for (int iIndex = 0; iIndex < arrObjects.Count; iIndex++)
			{
				Object3D obj3dCurrent = (Object3D)arrObjects[iIndex];
				if (obj3dCurrent == null)
					continue;
				// convert the object to eye coordinates
				double dtCurrent = obj3dCurrent.IntersectionTest(vectRayStart, vectRayEnd);
				// check for any intersection
				if (dtCurrent <= 0)
					continue;
				// intersection - ray has hit an object
				if (obj3dClosest == null || dtCurrent < dtClosest)
				{
					obj3dClosest = obj3dCurrent;
					dtClosest = dtCurrent;
				}
			}
			// no object intersect the ray, so return black, or 0 intensity
			if (obj3dClosest == null)
				return vectIntensity;
			// shadow ray flag
			bool bIsShadowRay = false;

			// establish a base intensity of the object
			vectIntensity = vectAmbientLightIntensity * vectAmbientLightColor * obj3dClosest.AmbientColor;

			// now, go through all point light sources, and see if it
			// influences the color of the light

			// get 3d position where the ray touched the object
			Vector vectRay = (vectRayEnd - vectRayStart).Normalize();
			Vector vectSurfacePoint = vectRayStart + dtClosest * vectRay;
			// get the normal of the surface hit
			Vector vectSurfaceNormal = obj3dClosest.GetSurfaceNormal(vectSurfacePoint);
			// tiny adjustment
			vectSurfacePoint += (vectSurfaceNormal * 0.01);
			for (int iIndex = 0; iIndex < arrPointLightSources.Count; iIndex++)
			{
				// get the point light source in the scene
				PointLightSource plsLight = (PointLightSource)arrPointLightSources[iIndex];
				if (plsLight == null)
					continue;
				bIsShadowRay = false;
				// convert the light source to UVW coodinates
				//Vector vectLightRay = plsLight.Position - vectSurfacePoint;
				// test to see how far away the light intersects the vector
				// dt = obj3dClosest.IntersectionTest(plsLight.Position, vectSurfacePoint);
				// test to see if the ray is in the shadow of another object wrt
				// the light source - assuming that the light ray to the surface of
				// of the object is 1.0, then any objects between the light source
				// and the current object will have a value between 0 and 1
				for (int iIndex2 = 0; iIndex2 < arrObjects.Count; iIndex2++)
				{
					Object3D obj3dCurrent = (Object3D)arrObjects[iIndex2];
					if (obj3dCurrent == null || obj3dClosest == obj3dCurrent)
						continue;
					dtClosest = obj3dCurrent.IntersectionTest(vectSurfacePoint, plsLight.Position);
					// yes, there is another object in the way of this object, so the
					// ray cast by the light source is actually a shadow ray for this object

					if (dtClosest > 0)
					{
						bIsShadowRay = true;
						break;
					}
				}
				// if the ray is a shadow ray, then the light from the light source does not contribute
				// to the light in the ray, so just continue on to the next light source
				if (bIsShadowRay)
					continue;

				// ok, so the light source contributes to the light shown on the object
				Vector vectLightVector = (plsLight.Position - vectSurfacePoint).Normalize();

				//Vector vectLightVector = new Vector(0, 1, 0);
				// calculate diffuse component of the light source
				Vector vectDiffuse = (plsLight.Diffuse * obj3dClosest.DiffuseColor *
					vectSurfaceNormal.DotProduct(vectLightVector)).Limit(0,1);
				vectIntensity += vectDiffuse;
				// get the inverse of the viewing ray

				// calculate specular component of the light source
				// get the "Halfway" vector
				Vector vectH = (vectLightVector + vectRay * -1).Normalize();
				vectSpecular = (plsLight.Specular * obj3dClosest.SpecularColor *
					Math.Pow(vectH.DotProduct(vectSurfaceNormal), obj3dClosest.SpecularPow)).Limit(0, 1);
				vectIntensity += vectSpecular;

			}
			Vector vectRecursiveRay = vectRay - (2 * vectRay.DotProduct(vectSurfaceNormal) * vectSurfaceNormal);
			Vector vectReflectiveIntensity = (TraceRay(vectSurfacePoint, vectRecursiveRay + vectSurfacePoint, 
				++iRayCount)).Limit(0, 1);
			vectIntensity += obj3dClosest.SpecularColor * vectReflectiveIntensity;

			// add the final and reflective colors together
			return vectIntensity; //  +vectReflectiveIntensity;
		}

		/// <summary>
		/// Adds a sphere to the ray tracer.
		/// </summary>
		/// <param name="dXPos"></param>
		/// <param name="dYPos"></param>
		/// <param name="dZPos"></param>
		/// <param name="clColor"></param>
		/// <param name="dRadius"></param>
		public Sphere AddSphere(double dXPos, double dYPos, double dZPos, double dRadius)
		{
			if (arrObjects == null)
				arrObjects = new ArrayList();
			Sphere sphNewSphere = new Sphere(dXPos, dYPos, dZPos, dRadius);
			arrObjects.Add(sphNewSphere);
			return sphNewSphere;
		}


		public void AddPointLightSource(double dXPos, double dYPos, double dZPos)
		{
			if (arrPointLightSources == null)
				arrPointLightSources = new ArrayList();
			arrPointLightSources.Add(new PointLightSource(dXPos, dYPos, dZPos));
		}


		private Vector ColorToVector(Color clColor){
			Vector vectColor = new Vector(3, true);
			vectColor[0] = clColor.R / 255.0;		// red component
			vectColor[1] = clColor.G / 255.0;		// green component
			vectColor[2] = clColor.B / 255.0;		// blue component
			return vectColor;
		}

		private Color VectorToColor(Vector vectColor)
		{
			vectColor = vectColor.Limit(0, 1.0);
			Color clNewColor = Color.FromArgb((int)(vectColor[0] * 255), (int)(vectColor[1] * 255.0),
								(int)(vectColor[2] * 255.0));
			return clNewColor;
		}

		
	}

	public class Camera
	{
		protected Vector vectPosition;		
		protected Vector vectLookAt;			
		protected Vector vectUp;			

		public Camera(double dXPos, double dYPos, double dZPos,
			double dLookAtX, double dLookAtY, double dLookAtZ,
			double dUpX, double dUpY, double dUpZ)
		{
			vectPosition = new Vector(dXPos, dYPos, dZPos);
			vectLookAt = new Vector(dLookAtX, dLookAtY, dLookAtZ);
			vectUp = new Vector(dUpX, dUpY, dUpZ);
		}

		public Vector Position
		{
			get
			{
				return vectPosition;
			}
		}

		/// <summary>
		/// Gets the U basis vector for the camera. (X)
		/// </summary>
		public Vector U
		{
			get
			{
				// the U basis vector of the camera is
				// the normalized cross product of the up 
				// and W vectors

				return MatrixMath.CrossProduct(vectUp, W);
			}
		}

		/// <summary>
		/// Gets the V basic vector of the camera (up)
		/// </summary>
		public Vector V
		{
			get
			{
				return MatrixMath.CrossProduct(W, U);
			}
		}

		/// <summary>
		/// Gets the W basis vector of the camera is a normalized
		/// gaze vector.
		/// </summary>
		public Vector W
		{
			get
			{
				return (vectLookAt - vectPosition).Normalize() * -1;
			}
		}
	}

	public abstract class Object3D
	{
		protected Vector vectPosition;

		protected Vector vectAmbientColor = new Vector(0.17, 0.1, 0.01);
		protected Vector vectDiffuseColor = new Vector(0.45, 0.29, 0.04);
		protected Vector vectSpecularColor = new Vector(0.9, 0.9, 0.2);
		protected double dSpecularPow = 35;
		protected double dReflection = 0.4;

		public Object3D(double dXPos, double dYPos, double dZPos)
		{
			vectPosition = new Vector(dXPos, dYPos, dZPos);
		}


		/// <summary>
		/// Resturns the X Y Z position of the object as a vector of size = 3.
		/// </summary>
		public Vector Position
		{
			get
			{
				return vectPosition;
			}
		}

		public Vector AmbientColor
		{
			get
			{
				return vectAmbientColor;
			}
			set
			{
				vectAmbientColor = new Vector(value);
			}
		}

		public Vector DiffuseColor
		{
			get
			{
				return vectDiffuseColor;
			}
			set
			{
				vectDiffuseColor = new Vector(value);
			}
		}

		public Vector SpecularColor
		{
			get
			{
				return vectSpecularColor;
			}
			set
			{
				vectSpecularColor = new Vector(value);
			}
		}

		public double SpecularPow
		{
			get
			{
				return dSpecularPow;
			}
			set
			{
				if (value < 0)
					value = 0;
				dSpecularPow = value;
			}
		}

		public double Reflection
		{
			get
			{
				return dReflection;
			}
			set
			{
				if (value < 0)
					value = 0;
				if (value > 1.0)
					value = 1.0;
				dReflection = value;
			}
		}

		/// <summary>
		/// An interstion test performs an analysis along a ray to determine if
		/// the ray passes through a 3d object in a scene. The value returned
		/// reflect the point along that ray where the insesction begins.
		/// 
		/// If zero is returned, then no intersection has occured.
		/// 
		/// To get the actual 3D position where the intersection occured, use the
		/// following formula (where t is a possible value returned):
		/// 
		/// p(t) = o + t(s - o)
		/// 
		/// Where o is the origin position (vectOrigin), t is the value returned 
		/// by this method, and s is the direction vector (vectRay).
		/// </summary>
		public abstract double IntersectionTest(Vector vectRayStart, Vector vectRayEnd);

		public abstract Vector GetSurfaceNormal(Vector vectSurfacePosition);
	}

	public class Sphere : Object3D
	{
		protected double dRadius;

		public Sphere(double dXPos, double dYPos, double dZPos, double dRadius): 
			base(dXPos, dYPos, dZPos)
		{
			this.dRadius = dRadius;
		}

		public double Radius
		{
			get
			{
				return dRadius;
			}
		}

		public override double IntersectionTest(Vector vectRayStart, Vector vectRayEnd)
		{
			Vector vectRayDirection = (vectRayEnd - vectRayStart).Normalize();
			double dA = vectRayDirection.DotProduct(vectRayDirection);
			Vector vectSphereToRayOrigin = vectRayStart - vectPosition;
			double dB = 2 * vectRayDirection.DotProduct(vectSphereToRayOrigin);
			double dC = (vectSphereToRayOrigin).DotProduct(vectSphereToRayOrigin) - (dRadius * dRadius);
			double dDiscriminant = (dB * dB) - (4 * dA * dC);

			double dt = 0;
			double dt1, dt2;
			if (dDiscriminant >= 0)  // has a solution
			{
				dt1 = (-dB - Math.Sqrt(dDiscriminant)) / (2.0 * dA);
				dt2 = (-dB + Math.Sqrt(dDiscriminant)) / (2.0 * dA);
				if (dt1 < dt2)
					dt = dt1;
				else
					dt = dt2;
			}
			return dt;
		}

		public override Vector GetSurfaceNormal(Vector vectSurfacePosition)
		{
			return (vectSurfacePosition - vectPosition).Normalize();
		}
	}

	public class Triangle : Object3D
	{
		Vector vectP1;
		Vector vectP2;
		Vector vectP3;
		double dSize;

		public Triangle(double dXPos, double dYPos, double dZPos, 
			double dSize): base(dXPos, dYPos, dZPos)
		{
			this.dSize = dSize;
		}

		public override double IntersectionTest(Vector vectRay, Vector vectOrigin)
		{
			throw new Exception("The method or operation is not implemented.");
		}

		public override Vector GetSurfaceNormal(Vector vectSurfacePosition)
		{
			Vector vectNormal = (vectSurfacePosition - vectPosition).Normalize(); // *(1.0 / dSize);
			return vectNormal;
		}
	}

	public class PointLightSource
	{
		private Vector vectPosition;

		private Vector vectDiffuse = new Vector(0.7, 0.7, 0.7);
		private Vector vectSpecular = new Vector(0.7, 0.7, 0.7);

		public PointLightSource(double dXPos, double dYPos, double dZPos)
		{
			vectPosition = new Vector(dXPos, dYPos, dZPos);
		}

		public Vector Position
		{
			get
			{
				return vectPosition;
			}
			set
			{
				vectPosition = new Vector(value);
			}
		}

		public Vector Diffuse
		{
			get
			{
				return vectDiffuse;
			}
		}

		public Vector Specular
		{
			get
			{
				return vectSpecular;
			}
		}
	}


}
