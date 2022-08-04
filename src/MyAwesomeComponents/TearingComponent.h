#pragma once
#include <MyAwesomeComponents/config.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/VecTypes.h>

#include <SofaMiscFem/TriangleFEMForceField.h>

#include <fstream>

namespace sofa::component::controller
{
	/**
	 * @brief TearingController Class
	 *
	 * Provides a Mouse & Keyboard user control on a Mechanical State.
	 * On a Rigid Particle, relative and absolute control is available.
	 */
	template<class DataTypes>
	class TearingComponent : public sofa::core::objectmodel::BaseObject
	{
	public:
		SOFA_CLASS(SOFA_TEMPLATE(TearingComponent, DataTypes), sofa::core::objectmodel::BaseObject);

		typedef defaulttype::Vec3Types DataTypes;
		typedef DataTypes::Coord Coord;
		typedef DataTypes::Real Real;

		/// Sofa API init method of the component
		void init() override;
		/// Sofa API reset method of the component
		void reset() override;

		/// Impl method that will compute the intersection and check if some element have to be removed.
		void doTear();

	protected:

		// Pointer to the target object collision model
		std::vector<core::CollisionModel*> m_surfaceCollisionModels;

		/// Default constructor
		TearingComponent();

		/// Default destructor
		~TearingComponent() override;


	public:
		Data<Real> d_yieldStress;
		sofa::component::forcefield::TriangleFEMForceField<DataTypes>* m_triangleFEMForceField;
		//sofa::component::forcefield::TriangularFEMForceField<DataTypes>* m_triangularFEMForceField;

		//SingleLink<TearingComponent, sofa::component::forcefield::TriangularFEMForceField<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_triangularForceField;
	};
}