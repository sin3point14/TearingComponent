#pragma once
#include <MyAwesomeComponents/config.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/VecTypes.h>

#include <SofaMiscFem/TriangularFEMForceField.h>
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>

#include <SofaUserInteraction/TopologicalChangeManager.h>

#include <SofaBoundaryCondition/FixedConstraint.h>

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
	class TearingComponent : public sofa::component::forcefield::TriangularFEMForceField<DataTypes>
	{
	public:
		SOFA_CLASS(SOFA_TEMPLATE(TearingComponent, DataTypes), SOFA_TEMPLATE(sofa::component::forcefield::TriangularFEMForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
	typedef typename sofa::defaulttype::Vec3Types Vec3types;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef sofa::core::topology::BaseMeshTopology::Index Index;
    typedef sofa::core::topology::BaseMeshTopology::Triangle Element;
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles VecElement;
    typedef sofa::core::topology::BaseMeshTopology::TrianglesAroundVertex TrianglesAroundVertex;

    typedef sofa::type::Quat<Real> Quat;
		// typedef defaulttype::Vec3Types DataTypes;
		// typedef typename DataTypes::Coord Coord;
		// typedef typename DataTypes::Real Real;

		/// Sofa API init method of the component
		void init() override;
		/// Sofa API reset method of the component
		// void reset() override;

		/// Impl method that will compute the intersection and check if some element have to be removed.
		// void doTear();
		Coord m_maxPrincipalStressDir;
		Real m_maxPrincipalStress;
		Size m_maxPrincipalStressIdx;

	protected:

		// Pointer to the target object collision model
		// std::vector<core::CollisionModel*> m_surfaceCollisionModels;
		sofa::component::topology::TriangleSetTopologyModifier* m_triangleMod;
		sofa::component::topology::TriangleSetTopologyContainer* m_triangleCon;
		sofa::component::projectiveconstraintset::FixedConstraint<sofa::defaulttype::Vec3Types>* m_fixedConstraint;

		sofa::component::topology::TriangleSetGeometryAlgorithms<sofa::defaulttype::Vec3Types>* m_triangleGeo;
		sofa::component::collision::TriangleCollisionModel<sofa::defaulttype::Vec3Types>* m_body;

		sofa::component::collision::TopologicalChangeManager topologyChangeManager;
		//sofa::component::collision::TriangleCollisionModel m_triangleCollisionModel;



		/// Default constructor
		// TearingComponent();

		/// Default destructor
		// ~TearingComponent() override;


	public:
    	void addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
		void draw(const core::visual::VisualParams* vparams) override;

		// Data<Real> d_yieldStress;
		// sofa::component::forcefield::TriangularFEMForceField<DataTypes>* m_triangleFEMForceField;
		//sofa::component::forcefield::TriangularFEMForceField<DataTypes>* m_triangularFEMForceField;

		//SingleLink<TearingComponent, sofa::component::forcefield::TriangularFEMForceField<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_triangularForceField;
	};
}