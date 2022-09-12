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
	core::topology::BaseMeshTopology::TriangleID TearingComponent<DataTypes>::findCutEndpoint(core::topology::BaseMeshTopology::TriangleID source,
		Coord maxPrincipalStressDir,
		const VecCoord& points,
		bool alongPosX)
	{
		core::topology::BaseMeshTopology::TriangleID aIndex = source;

		Coord aPoints[3];
		sofa::type::Vec3 aCentroid;

		m_triangleGeo->getTriangleVertexCoordinates(aIndex, aPoints);

		aCentroid.clear();
		aCentroid = (aPoints[0] + aPoints[1] + aPoints[2]) / 3.0;

		sofa::type::vector< sofa::core::topology::TopologyElementType> topoPathList;
		sofa::type::vector<Index> indicesList;
		sofa::type::vector< sofa::type::Vec<3, double> > coords2List;

		core::topology::BaseMeshTopology::EdgesInTriangle edges = m_triangleCon->getEdgesInTriangle(aIndex);

		const float planeA = maxPrincipalStressDir[0];
		const float planeB = maxPrincipalStressDir[1];
		const float planeC = maxPrincipalStressDir[2];
		const float planeD = -planeA * aCentroid[0] - planeB * aCentroid[1] - planeC * aCentroid[2];
		const auto planeEqn = [&](Coord p) {
			return planeA * p[0] + planeB * p[1] + planeC * p[2] + planeD;
		};

		core::topology::BaseMeshTopology::TriangleID endPointID = -1;
		// distance of centroid from source
		float endPointDistance = 0;
		while (endPointDistance < 40.0f) {
			core::topology::BaseMeshTopology::TriangleID prevNext = endPointID;
			for (auto&& currID : m_triangleCon->getElementAroundElement(aIndex)) {
				core::topology::BaseMeshTopology::PointID toBeCut1, toBeCut2;
				auto& t = m_triangleCon->getTriangle(currID);
				Coord p1 = points[t[0]];
				Coord p2 = points[t[1]];
				Coord p3 = points[t[2]];
				float dist1 = (aCentroid - p1).norm2();
				float dist2 = (aCentroid - p2).norm2();
				float dist3 = (aCentroid - p3).norm2();

				bool pos1 = planeEqn(p1) > 0.0f;
				bool pos2 = planeEqn(p2) > 0.0f;
				bool pos3 = planeEqn(p3) > 0.0f;
				if ((pos1 && pos2 && pos3) || (!pos1 && !pos2 && !pos3))
					continue;
				if (pos1 == pos2) {
					toBeCut1 = t[2];
					if (dist1 > dist2) {
						toBeCut2 = t[0];
					}
					else {
						toBeCut2 = t[1];
					}
				}
				else if (pos2 == pos3) {
					toBeCut1 = t[0];
					if (dist2 > dist3) {
						toBeCut2 = t[1];
					}
					else {
						toBeCut2 = t[2];
					}
				}
				else {
					toBeCut1 = t[1];
					if (dist3 > dist1) {
						toBeCut2 = t[2];
					}
					else {
						toBeCut2 = t[0];
					}
				}
				core::topology::BaseMeshTopology::EdgeID cutEdge = m_triangleCon->getEdgeIndex(toBeCut1, toBeCut2);
				core::topology::BaseMeshTopology::TrianglesAroundEdge tris = m_triangleCon->getTrianglesAroundEdge(cutEdge);
				auto other = tris[0] == currID ? tris[1] : tris[0];
				Coord otherCoord[3];
				m_triangleGeo->getTriangleVertexCoordinates(other, otherCoord);
				sofa::type::Vec3 otherCentroid = (otherCoord[0] + otherCoord[1] + otherCoord[2]) / 3.0;
				float distCentroid = (aCentroid - otherCentroid).norm();
				float xDiffPos = (aCentroid - otherCentroid).x() > 0;
				if (distCentroid > endPointDistance && xDiffPos == alongPosX) {
					endPointDistance = distCentroid;
					endPointID = other;
				}
			}
			if (prevNext == endPointID)
				break;
		}
		std::cout << endPointDistance << std::endl;
		return endPointID;
	}

	template <class DataTypes>
	void TearingComponent<DataTypes>::addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) 
	{
		// support only one cut as of now
		static bool onlyOnce = true;

		Size nbTriangles = this->m_topology->getNbTriangles();
    	const auto& triangleInf = this->triangleInfo.getValue();

		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::addForce(mparams, f, x, v);

		const VecCoord& x1 = this->mstate->read(core::ConstVecCoordId::position())->getValue();
		m_maxPrincipalStress = 0;

        for(Size i=0; i< nbTriangles; ++i)
        {
			Deriv direction;
			Real value;
			// doesn't actually do any fracture stuff
			sofa::component::forcefield::TriangularFEMForceField<DataTypes>::getFractureCriteria(i, direction, value);
			direction.normalize();
			if (value > m_maxPrincipalStress)
			{
				m_maxPrincipalStress = value;
				m_maxPrincipalStressDir = direction;
				m_maxPrincipalStressIdx = i;
			}
			
        }

		// heuristic
		if(m_maxPrincipalStress > 2.0f && onlyOnce)
		{
			onlyOnce = false;
			// PointIDs -> Points
			const typename DataTypes::VecCoord& points = m_triangleGeo->getDOF()->read(core::ConstVecCoordId::position())->getValue();

			core::topology::BaseMeshTopology::TriangleID tri1 = findCutEndpoint(m_maxPrincipalStressIdx, m_maxPrincipalStressDir, points, true);
			core::topology::BaseMeshTopology::TriangleID tri2 = findCutEndpoint(m_maxPrincipalStressIdx, -m_maxPrincipalStressDir, points, false);

			if (tri1 != -1 && tri2 != -1) {
				sofa::type::Vec3 tri1Centroid;
				Coord tri1Points[3];
				m_triangleGeo->getTriangleVertexCoordinates(tri1, tri1Points);

				tri1Centroid = (tri1Points[0] + tri1Points[1] + tri1Points[2]) / 3.0;

				sofa::type::Vec3 tri2Centroid;
				Coord tri2Points[3];
				m_triangleGeo->getTriangleVertexCoordinates(tri2, tri2Points);

				tri2Centroid = (tri2Points[0] + tri2Points[1] + tri2Points[2]) / 3.0;

				m_topologyChangeManager.incisionCollisionModel(m_collisionModel, tri1, tri1Centroid, m_collisionModel, tri2, tri2Centroid, 50);
			}
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

	template <class DataTypes>
	void TearingComponent<DataTypes>::init()
	{
		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::init();
		this->m_topology->getContext()->get(m_triangleCon);
		this->m_topology->getContext()->get(m_triangleGeo);
		this->m_topology->getContext()->get(m_fixedConstraint);
		this->m_topology->getContext()->get(m_collisionModel);

		sofa::type::vector<core::topology::BaseMeshTopology::PointID> borderPoints = m_triangleCon->getPointsOnBorder();

		// fixing all border points to mimic a membrane
		for (auto&& pid : borderPoints) {
			m_fixedConstraint->addConstraint(pid);
		}

		if (m_triangleGeo == nullptr || m_collisionModel == nullptr || m_triangleCon == nullptr)
		{
			msg_error() << "No TriangleSetTopologyContainer, FixedConstraint, TriangleSetGeometryAlgorithms, TopologicalChangeManager or TriangleCollisionModel not found";
			sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
			return;
		}
	}
}
