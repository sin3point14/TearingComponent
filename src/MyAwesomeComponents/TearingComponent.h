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

		/// Sofa API init method of the component
		void init() override;

		Coord m_maxPrincipalStressDir;
		Real m_maxPrincipalStress;
		Size m_maxPrincipalStressIdx;

	protected:

		class TriangleData {
		public:
			core::topology::BaseMeshTopology::TriangleID t;
			core::topology::BaseMeshTopology::PointID p1, p2, p3;
			core::topology::BaseMeshTopology::EdgeID e12, e23, e31;
			TriangleData(core::topology::BaseMeshTopology::TriangleID tID,
				const VecCoord& points,
				sofa::component::topology::TriangleSetTopologyContainer* triangleCon)
				: t(tID)
			{
				sofa::topology::Triangle t = triangleCon->getTriangle(tID);
				p1 = t[0];
				p2 = t[1];
				p3 = t[2];
				e12 = triangleCon->getEdgeIndex(p1, p2);
				e23 = triangleCon->getEdgeIndex(p2, p3);
				e31 = triangleCon->getEdgeIndex(p3, p1);
			}
		};

		sofa::component::topology::TriangleSetTopologyContainer* m_triangleCon;
		sofa::component::projectiveconstraintset::FixedConstraint<sofa::defaulttype::Vec3Types>* m_fixedConstraint;
		sofa::component::topology::TriangleSetGeometryAlgorithms<sofa::defaulttype::Vec3Types>* m_triangleGeo;
		sofa::component::topology::TriangleSetTopologyModifier* m_triangleMod;
		sofa::component::collision::TriangleCollisionModel<sofa::defaulttype::Vec3Types>* m_collisionModel;
		sofa::component::collision::TopologicalChangeManager m_topologyChangeManager;

		core::topology::BaseMeshTopology::TriangleID getOtherTriangle(const TriangleData& curr,
			core::topology::BaseMeshTopology::EdgeID edge);

		// Adjacent Triangles sharing en edge
		sofa::type::vector<core::topology::BaseMeshTopology::TriangleID>& getAdjacentTriangles(const TriangleData& curr);

		float fixEdgeBarycentricParameter(core::topology::BaseMeshTopology::EdgeID edge,
			core::topology::BaseMeshTopology::PointID p1, float t);

		// finds one end point, need to call this twice to
		// with true and false along_posX values to get both
		// end points
		TriangleData findCutEndpoint(core::topology::BaseMeshTopology::TriangleID source,
			Coord maxPrincipalStressDir,
			const VecCoord& points,
			bool alongPosX,
			sofa::type::vector< sofa::core::topology::TopologyElementType>& topoPath_list,
			sofa::type::vector<Index>& indices_list,
			sofa::type::vector<Coord>& coords2_list);

	public:
    	void addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
		void draw(const core::visual::VisualParams* vparams) override;
	};
}