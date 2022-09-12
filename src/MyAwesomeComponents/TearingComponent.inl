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
	core::topology::BaseMeshTopology::TriangleID TearingComponent<DataTypes>::getOtherTriangle(const TriangleData& curr,
		core::topology::BaseMeshTopology::EdgeID edge)
	{
		core::topology::BaseMeshTopology::TrianglesAroundEdge tris = m_triangleCon->getTrianglesAroundEdge(edge);
		core::topology::BaseMeshTopology::TriangleID otherTri = tris[0] == curr.t ? tris[1] : tris[0];
		return otherTri;
	}

	template <class DataTypes>
	sofa::type::vector<core::topology::BaseMeshTopology::TriangleID>& TearingComponent<DataTypes>::getAdjacentTriangles(const TriangleData& curr)
	{
		sofa::type::vector<core::topology::BaseMeshTopology::TriangleID> toRet;
		toRet.push_back(getOtherTriangle(curr, curr.e12));
		toRet.push_back(getOtherTriangle(curr, curr.e23));
		toRet.push_back(getOtherTriangle(curr, curr.e31));
		return toRet;
	}

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
		sofa::type::vector<sofa::type::Vec<3, double>> coords2List;

		core::topology::BaseMeshTopology::EdgesInTriangle edges = m_triangleCon->getEdgesInTriangle(aIndex);

		const float planeA = maxPrincipalStressDir[0];
		const float planeB = maxPrincipalStressDir[1];
		const float planeC = maxPrincipalStressDir[2];
		const float planeD = -planeA * aCentroid[0] - planeB * aCentroid[1] - planeC * aCentroid[2];
		// returns t for intersection of plane and p1 + (p2 - p1) * t
		const auto planeLineIntersction = [&](Coord p1, Coord p2) {
			return (-planeA * p1[0] - planeB * p1[1] - planeC * p1[2] - planeD) 
				/ (planeA * (p2[0] - p1[0]) + planeB * (p2[1] - p1[1]) + planeC * (p2[2] - p1[2]));
		};

		// stores all triangles that will be cut
		sofa::type::vector<core::topology::BaseMeshTopology::TriangleID> incisionIDs = { source };
		
		// distance of centroid from source
		float nextCandidateDistance = 0;
		core::topology::BaseMeshTopology::TriangleID nextCandidate;

		while (nextCandidateDistance < 40.0f) {
			nextCandidate = sofa::InvalidID;
			TriangleData lastCutTri(incisionIDs[incisionIDs.size() - 1], points, m_triangleCon);

			for (auto&& currID : getAdjacentTriangles(lastCutTri)) {
				core::topology::BaseMeshTopology::EdgeID toBeCut;
				TriangleData currTri(currID, points, m_triangleCon);
				Coord p1 = points[currTri.p1];
				Coord p2 = points[currTri.p2];
				Coord p3 = points[currTri.p3];
				float dist1 = (aCentroid - p1).norm2();
				float dist2 = (aCentroid - p2).norm2();
				float dist3 = (aCentroid - p3).norm2();

				float t12 = planeLineIntersction(p1, p2);
				float t23 = planeLineIntersction(p2, p3);
				float t31 = planeLineIntersction(p3, p1);

				bool t12Intersection = t12 > 0.0f && t12 < 1.0f;
				bool t23Intersection = t23 > 0.0f && t23 < 1.0f;
				bool t31Intersection = t31 > 0.0f && t31 < 1.0f;

				if (!(t12Intersection || t23Intersection || t31Intersection))
					continue;

				// determine which edge to cur
				// first we find the edge between current and previous endPoint
				// as there is no need to consider it
				// take an assumption that there are no looping cuts i.e. 
				// out of the adjacent triangles only one is present inincisionIDs
				// which will always be the last element of incisionIDs
				core::topology::BaseMeshTopology::TriangleID otherTri12 = getOtherTriangle(currTri, currTri.e12);
				core::topology::BaseMeshTopology::TriangleID otherTri23 = getOtherTriangle(currTri, currTri.e23);
				core::topology::BaseMeshTopology::TriangleID otherTri31 = getOtherTriangle(currTri, currTri.e31);

				core::topology::BaseMeshTopology::TriangleID last = incisionIDs[incisionIDs.size() - 1];

				if (otherTri12 == lastCutTri.t) {
					if (t23Intersection) {
						toBeCut = currTri.e23;
					}
					else if (t31Intersection) {
						toBeCut = currTri.e31;
					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else if (otherTri23 == lastCutTri.t) {
					if (t31Intersection) {
						toBeCut = currTri.e31;
					}
					else if (t12Intersection) {
						toBeCut = currTri.e12;
					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else if (otherTri31 == lastCutTri.t) {
					if (t12Intersection) {
						toBeCut = currTri.e12;
					}
					else if (t23Intersection) {
						toBeCut = currTri.e23;
					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else {
					throw std::exception("Unreachable!");
				}
				
				core::topology::BaseMeshTopology::TriangleID other = getOtherTriangle(currTri, toBeCut);
				Coord otherCoord[3];
				m_triangleGeo->getTriangleVertexCoordinates(other, otherCoord);
				sofa::type::Vec3 otherCentroid = (otherCoord[0] + otherCoord[1] + otherCoord[2]) / 3.0;
				float distCentroid = (aCentroid - otherCentroid).norm();
				float xDiffPos = (aCentroid - otherCentroid).x() > 0;
				if (distCentroid > nextCandidateDistance && xDiffPos == alongPosX) {
					nextCandidateDistance = distCentroid;
					nextCandidate = other;
				}
			}
			if (nextCandidate == sofa::InvalidID)
				break;
			incisionIDs.push_back(nextCandidate);
		}
		std::cout << nextCandidateDistance << std::endl;
		return incisionIDs[incisionIDs.size() - 1];;
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
