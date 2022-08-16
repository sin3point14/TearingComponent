#include "TearingComponent.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.h>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/ColorMap.h>
#include <sofa/type/RGBAColor.h>

#include <sofa/core/topology/TopologyData.inl>

#include <Eigen/SVD>
#include <limits>

namespace sofa::component::controller
{
	template <class DataTypes>
	void TearingComponent<DataTypes>::addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) 
	{
		static bool onlyOnce = true;
		Size nbTriangles = this->m_topology->getNbTriangles();

    	const auto& triangleInf = this->triangleInfo.getValue();

		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::addForce(mparams, f, x, v);
		const VecCoord& x1 = this->mstate->read(core::ConstVecCoordId::position())->getValue();
		m_maxPrincipalStress = 0;
        // std::vector<sofa::type::Vector3> vertices;
        for(Size i=0; i< nbTriangles; ++i)
        {
			Deriv direction;
			Real value;
			sofa::component::forcefield::TriangularFEMForceField<DataTypes>::getFractureCriteria(i, direction, value);
			direction.normalize();
			if (value > m_maxPrincipalStress)
			{
				m_maxPrincipalStress = value;
				m_maxPrincipalStressDir = direction;
				m_maxPrincipalStressIdx = i;
			}
			
            // const sofa::core::topology::BaseMeshTopology::Triangle& tri = this->m_topology->getTriangle(i);
            // Index a = tri[0];
            // Index b = tri[1];
            // Index c = tri[2];
            // Coord center = (x1[a]+x1[b]+x1[c])/3;
            // Coord d = triangleInf[i].principalStressDirection*2.5; //was 0.25
            // vertices.push_back(sofa::type::Vector3(center));
            // vertices.push_back(sofa::type::Vector3(center+d));
        }
        // vparams->drawTool()->drawLines(vertices, 1, sofa::type::RGBAColor(1, 0, 1, 1));

		if(m_maxPrincipalStress > 100.0f && onlyOnce)
		{
   // const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

			//const sofa::core::topology::BaseMeshTopology::Triangle& tri = this->m_topology->getTriangle(m_maxPrincipalStressIdx);
			//Index a = tri[0];
			//Index b = tri[1];
			//Index c = tri[2];
			//Coord normal = cross(x[a] - x[b], x[a] - x[c]);
   //         Coord center = (x[a]+x[b]+x[c])/3;
   //         Coord d = m_maxPrincipalStressDir*2.5; //was 0.25
			//Coord tearingDirection = cross(normal, m_maxPrincipalStressDir);
   //         Coord d2 = tearingDirection*2.5; //was 0.25

			//Coord p1, p2;

			//Coord s1 = (x[a] - x[b]);
			//Coord s2 = (x[b] - x[c]);
			//Coord s3 = (x[c] - x[a]);
			//s1.normalize();
			//s2.normalize();
			//s3.normalize();

			//double dot1 = fabs(dot(s1, tearingDirection));
			//double dot2 = fabs(dot(s2, tearingDirection));
			//double dot3 = fabs(dot(s3, tearingDirection));

			//if(dot1 <= dot2)
			//{
			//	if(dot1 <= dot3)
			//	{
			//		p1 = (x[b] + x[c]) / 2;
			//		p2 = (x[c] + x[a]) / 2;
			//	}
			//	else
			//	{
			//		p1 = (x[a] + x[b]) / 2;
			//		p2 = (x[b] + x[c]) / 2;
			//	}
			//}
			//else
			//{
			//	if(dot2 <= dot3)
			//	{
			//		p1 = (x[a] + x[b]) / 2;
			//		p2 = (x[c] + x[a]) / 2;
			//	}
			//	else
			//	{
			//		p1 = (x[a] + x[b]) / 2;
			//		p2 = (x[b] + x[c]) / 2;
			//	}
			//}

			//sofa::type::vector<sofa::core::topology::TopologyElementType>       topoPath_list;
			//sofa::type::vector<Index> indices_list;
			//sofa::type::vector<Coord> coords2_list;

			//bool isPathOk =
			//	m_triangleGeo->computeIntersectedObjectsList(
			//			sofa::InvalidID, p1, p2, m_maxPrincipalStressIdx, m_maxPrincipalStressIdx,
			//			topoPath_list, indices_list,
			//			coords2_list);
			//if (!isPathOk)
   //             {
   //                 msg_error() << "Invalid path in computeIntersectedPointsList";
   //                 return;
   //             }
			//sofa::type::vector<Index> new_edges;
			//
			//m_triangleGeo->SplitAlongPath(sofa::InvalidID, p1, sofa::InvalidID, p2,
   //                     topoPath_list, indices_list, coords2_list,
   //                     new_edges, 0.1, 0.25);
			//sofa::type::vector<Index> new_points;
			//sofa::type::vector<Index> end_points;
			//bool reachBorder = false;
			//m_triangleGeo->InciseAlongEdgeList(new_edges,
   //                     new_points, end_points, reachBorder);
			onlyOnce = false;
			core::topology::BaseMeshTopology::TriangleID ind_a = 1670;
			core::topology::BaseMeshTopology::TriangleID ind_b = 1137;

			Coord aCoord[3];
			Coord bCoord[3];
			sofa::type::Vec3 a;
			sofa::type::Vec3 b;

			sofa::Index a_last = sofa::InvalidID;
			sofa::Index b_last = sofa::InvalidID;

			sofa::component::topology::TriangleSetTopologyModifier* triangleMod;
			m_topology->getContext()->get(triangleMod);

			sofa::component::topology::TriangleSetGeometryAlgorithms<sofa::defaulttype::Vec3Types>* triangleGeo;
			m_topology->getContext()->get(triangleGeo);

			triangleGeo->getTriangleVertexCoordinates(ind_a, aCoord);
			a.clear();
			a = (aCoord[0] + aCoord[1] + aCoord[2]) / 3.0;

			triangleGeo->getTriangleVertexCoordinates(ind_b, bCoord);
			b.clear();
			b = (bCoord[0] + bCoord[1] + bCoord[2]) / 3.0;

			sofa::type::vector< sofa::core::topology::TopologyElementType> topoPath_list;
			sofa::type::vector<Index> indices_list;
			sofa::type::vector< sofa::type::Vec<3, double> > coords2_list;

			//errorTrianglesIndices.push_back(ind_ta);
			//errorTrianglesIndices.push_back(ind_tb);

			topologyChangeManager.incisionCollisionModel(m_body, ind_a, a, m_body, ind_b, b, 50);

			//bool isPathOk = triangleGeo->computeIntersectedObjectsList(sofa::InvalidID, a, b, ind_a, ind_b, topoPath_list, indices_list, coords2_list);

			//if (!isPathOk)
			//{
			//	msg_error() << "While computing computeIntersectedPointsList between triangles '"
			//		 << "' and '" << "' at time = '" << getContext()->getTime() << "'";

			//	msg_error() << " a = " << a << " b = " << b << msgendl
			//		<< "ind_ta = " << ind_a << " ind_tb = " << ind_b;

			//	return;
			//}


			//sofa::type::vector< Index > new_edges;


			//SReal  m_epsilonSnapPath = 0.1; ///< epsilon snap path
			//SReal  m_epsilonSnapBorder = 0.25; ///< epsilon snap path

			////Split triangles to create edges along a path given as a the list of existing edges and triangles crossed by it.
			//triangleGeo->SplitAlongPath(a_last, a, b_last, b, topoPath_list, indices_list, coords2_list, new_edges, m_epsilonSnapPath, m_epsilonSnapBorder);

			//sofa::type::vector<Index> new_points;
			//sofa::type::vector<Index> end_points;
			//bool reachBorder = false;

			////Duplicates the given edges
			//triangleGeo->InciseAlongEdgeList(new_edges, new_points, end_points, reachBorder);

			//msg_info_when(reachBorder) << "Incision has reached a border.";

			//triangleMod->notifyEndingEvent();
		}
	}

	template<class DataTypes>
	void TearingComponent<DataTypes>::draw(const core::visual::VisualParams* vparams)
	{
		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::draw(vparams);
		if(m_maxPrincipalStress > 100.0f)
		{
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

			const sofa::core::topology::BaseMeshTopology::Triangle& tri = this->m_topology->getTriangle(m_maxPrincipalStressIdx);
			Index a = tri[0];
			Index b = tri[1];
			Index c = tri[2];
			Coord normal = cross(x[a] - x[b], x[a] - x[c]);
            Coord center = (x[a]+x[b]+x[c])/3;
            Coord d = m_maxPrincipalStressDir*2.5; //was 0.25
			Coord tearingDirection = cross(normal, m_maxPrincipalStressDir);
            Coord d2 = tearingDirection*2.5; //was 0.25
        	vparams->drawTool()->drawLine(center, center+d, sofa::type::RGBAColor(0, 0, 1, 1));
        	vparams->drawTool()->drawLine(center, center+d2, sofa::type::RGBAColor(1, 0, 0, 1));
		}
	}

	// template <class DataTypes>
	// TearingComponent<DataTypes>::TearingComponent()
	// 	: d_yieldStress(initData(&d_yieldStress, (Real)10.0f, "yieldStress", "Yield Stress"))
	// 	, m_triangleFEMForceField(nullptr)
	// {
	// 	this->f_listening.setValue(true);
	// }

	// template <class DataTypes>
	// TearingComponent<DataTypes>::~TearingComponent()
	// {
	// }


	template <class DataTypes>
	void TearingComponent<DataTypes>::init()
	{
		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::init();
		// We look for a CollisionModel identified with the CarvingSurface Tag.

		// m_triangleFEMForceField = getContext()->get<sofa::component::forcefield::TriangularFEMForceField<DataTypes>>();
		// sofa::component::topology::TriangleSetTopologyModifier* triangleMod;
		this->m_topology->getContext()->get(m_triangleMod);

		// sofa::component::topology::TriangleSetGeometryAlgorithms<Vec3Types>* triangleGeo;
		this->m_topology->getContext()->get(m_triangleGeo);
		this->m_topology->getContext()->get(m_body);


		if (m_triangleMod == nullptr || m_triangleGeo == nullptr || m_body == nullptr)
		{
			msg_error() << "No TriangleSetTopologyModifier, TriangleSetGeometryAlgorithms or Collision model not found";
			sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
			return;
		}

		// if (m_triangleFEMForceField == nullptr)
		// {
		// 	msg_error() << "No TriangularFEMForceField component found";
		// 	sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
		// 	return;
		// }
	}


	// template <class DataTypes>
	// void TearingComponent<DataTypes>::reset()
	// {
	// }


	// template <class DataTypes>
	// void TearingComponent<DataTypes>::doTear()
	// {
	// 	if (m_triangleFEMForceField)
	// 	{

	// 	}
	// }

}
