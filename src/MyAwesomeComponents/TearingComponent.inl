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
	float TearingComponent<DataTypes>::fixEdgeBarycentricParameter(core::topology::BaseMeshTopology::EdgeID edge,
		core::topology::BaseMeshTopology::PointID p1, float t)
	{
		if (m_triangleCon->getEdge(edge)[0] == p1)
		{
			return t;
		}
		else
		{
			return 1 - t;
		}
	}

	template <class DataTypes>
	typename TearingComponent<DataTypes>::TriangleData TearingComponent<DataTypes>::findCutEndpoint(core::topology::BaseMeshTopology::TriangleID source,
		Coord maxPrincipalStressDir,
		const VecCoord& points,
		bool alongPosX,
		sofa::type::vector< sofa::core::topology::TopologyElementType>& topoPath_list,
		sofa::type::vector<Index>& indices_list,
		sofa::type::vector<Coord>& coords2_list)
	{
		core::topology::BaseMeshTopology::TriangleID aIndex = source;

		Coord aPoints[3];
		sofa::type::Vec3 aCentroid;

		m_triangleGeo->getTriangleVertexCoordinates(aIndex, aPoints);

		aCentroid.clear();
		aCentroid = (aPoints[0] + aPoints[1] + aPoints[2]) / 3.0;

		sofa::type::vector< sofa::core::topology::TopologyElementType> topoPathList;
		sofa::type::vector<Index> indicesList;
		sofa::type::vector<Coord> coords2List;

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
			core::topology::BaseMeshTopology::EdgeID cutEdgeFar = sofa::InvalidID;
			core::topology::BaseMeshTopology::EdgeID cutEdgeClose = sofa::InvalidID;
			float cutRatioFar = -1;
			float cutRatioClose = -1;

			for (auto&& currID : getAdjacentTriangles(lastCutTri)) {
				// TODO: process 1 triangle every iteration, this code has 
				// become extremly convoluted because I didn't realise I
				// was processing 2 triangles at once
				// This can be simplifying by just finding an edge that 
				// intersects the plane and selecting the other triangle
				// on the edge as the next triangle
				// 
				// Since my algortihm advances 2 triangles every iteration,
				// We calculate 2 edge cuts, this holds the farther one
				core::topology::BaseMeshTopology::EdgeID toBeCutFar;
				// This holds the near edge
				core::topology::BaseMeshTopology::EdgeID toBeCutClose;
				float toBeCutFarRatio = -1;
				float toBeCutCloseRatio = -1;
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

				// determine which edge to cut
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
					toBeCutClose = currTri.e12;
					toBeCutCloseRatio = fixEdgeBarycentricParameter(currTri.e12, currTri.p1, t12);
					if (t23Intersection) {
						toBeCutFar = currTri.e23;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e23, currTri.p2, t23);
					}
					else if (t31Intersection) {
						toBeCutFar = currTri.e31;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e31, currTri.p3, t31);

					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else if (otherTri23 == lastCutTri.t) {
					toBeCutClose = currTri.e23;
					toBeCutCloseRatio = fixEdgeBarycentricParameter(currTri.e23, currTri.p2, t23);
					if (t31Intersection) {
						toBeCutFar = currTri.e31;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e31, currTri.p3, t31);
					}
					else if (t12Intersection) {
						toBeCutFar = currTri.e12;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e12, currTri.p1, t12);
					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else if (otherTri31 == lastCutTri.t) {
					toBeCutClose = currTri.e31;
					toBeCutCloseRatio = fixEdgeBarycentricParameter(currTri.e31, currTri.p3, t31);
					if (t12Intersection) {
						toBeCutFar = currTri.e12;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e12, currTri.p1, t12);
					}
					else if (t23Intersection) {
						toBeCutFar = currTri.e23;
						toBeCutFarRatio = fixEdgeBarycentricParameter(currTri.e23, currTri.p2, t23);
					}
					else {
						throw std::exception("Unreachable!");
					}
				}
				else {
					throw std::exception("Unreachable!");
				}
				
				core::topology::BaseMeshTopology::TriangleID other = getOtherTriangle(currTri, toBeCutFar);
				Coord otherCoord[3];
				m_triangleGeo->getTriangleVertexCoordinates(other, otherCoord);
				sofa::type::Vec3 otherCentroid = (otherCoord[0] + otherCoord[1] + otherCoord[2]) / 3.0;
				float distCentroid = (aCentroid - otherCentroid).norm();
				float xDiffPos = (aCentroid - otherCentroid).x() > 0;
				if (distCentroid > nextCandidateDistance && xDiffPos == alongPosX) {
					nextCandidateDistance = distCentroid;
					nextCandidate = other;
					cutEdgeFar = toBeCutFar;
					cutEdgeClose = toBeCutClose;
					cutRatioFar = toBeCutFarRatio;
					cutRatioClose = toBeCutCloseRatio;
				}
			}
			if (nextCandidate == sofa::InvalidID)
				break;
			incisionIDs.push_back(nextCandidate);

			topoPath_list.push_back(core::topology::TopologyElementType::EDGE);
			indices_list.push_back(cutEdgeClose);
			topoPath_list.push_back(core::topology::TopologyElementType::EDGE);
			indices_list.push_back(cutEdgeFar);
			Coord baryCoords;
			baryCoords[0] = cutRatioClose;
			baryCoords[1] = 0.0;
			baryCoords[2] = 0.0;
			coords2_list.push_back(baryCoords);
			baryCoords[0] = cutRatioFar;
			coords2_list.push_back(baryCoords);
		}

		// End point: centroid of last triangle
		core::topology::BaseMeshTopology::TriangleID lastTri = incisionIDs[incisionIDs.size() - 1];
		Coord baryCoords;
		baryCoords[0] = 1./3;
		baryCoords[1] = 1./3;
		baryCoords[2] = 1./3;
		topoPath_list.push_back(core::topology::TopologyElementType::TRIANGLE);
		indices_list.push_back(lastTri);
		coords2_list.push_back(baryCoords);

		std::cout << nextCandidateDistance << std::endl;
		return TriangleData(incisionIDs[incisionIDs.size() - 1], points, m_triangleCon);
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

			sofa::type::vector< sofa::core::topology::TopologyElementType> topoPath_list;
			sofa::type::vector<Index> indices_list;
			sofa::type::vector<Coord> coords2_list;

			TriangleData tri1 = findCutEndpoint(m_maxPrincipalStressIdx, m_maxPrincipalStressDir, points, true, topoPath_list, indices_list, coords2_list);
			
			sofa::type::vector< sofa::core::topology::TopologyElementType> topoPath_list2;
			sofa::type::vector<Index> indices_list2;
			sofa::type::vector<Coord> coords2_list2;
			
			TriangleData tri2 = findCutEndpoint(m_maxPrincipalStressIdx, -m_maxPrincipalStressDir, points, false, topoPath_list2, indices_list2, coords2_list2);

			std::reverse(topoPath_list.begin(), topoPath_list.end());
			std::reverse(indices_list.begin(), indices_list.end());
			std::reverse(coords2_list.begin(), coords2_list.end());

			topoPath_list.insert(topoPath_list.end(), topoPath_list2.begin(), topoPath_list2.end());
			indices_list.insert(indices_list.end(), indices_list2.begin(), indices_list2.end());
			coords2_list.insert(coords2_list.end(), coords2_list2.begin(), coords2_list2.end());

			sofa::type::vector<core::topology::BaseMeshTopology::Edge> edges;

			for (auto&& e : indices_list)
			{
				std::cout << m_triangleCon->getEdge(e) <<  "  ";
			}

			std::cout << edges;

			sofa::type::vector< Index > new_edges;
			int result = m_triangleGeo->SplitAlongPath(sofa::InvalidID, (points[tri1.p1] + points[tri1.p2] + points[tri1.p3]) / 3,
				sofa::InvalidID, (points[tri2.p1] + points[tri2.p2] + points[tri2.p3]) / 3,
				topoPath_list, indices_list, coords2_list, new_edges, 0.0, 0.0);

			if (result == -1)
			{
				std::cout << "SplitAlongPath FAILED" << std::endl;
				return;
			}

			sofa::type::vector<Index> new_points;
			sofa::type::vector<Index> end_points;
			bool reachBorder = false;
			bool incision_ok = m_triangleGeo->InciseAlongEdgeList(new_edges, new_points, end_points, reachBorder);

			if (!incision_ok)
			{
				std::cout << "InciseAlongEdgeList FAILED" << std::endl;
				return;
			}

			m_triangleMod->notifyEndingEvent();
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
		this->m_topology->getContext()->get(m_triangleMod);
		this->m_topology->getContext()->get(m_fixedConstraint);
		this->m_topology->getContext()->get(m_collisionModel);

		sofa::type::vector<core::topology::BaseMeshTopology::PointID> borderPoints = m_triangleCon->getPointsOnBorder();

		// fixing all border points to mimic a membrane
		for (auto&& pid : borderPoints) {
			m_fixedConstraint->addConstraint(pid);
		}

		if (m_triangleGeo == nullptr || m_collisionModel == nullptr || m_triangleCon == nullptr)
		{
			msg_error() << "No TriangleSetTopologyContainer, FixedConstraint, TriangleSetGeometryAlgorithms, TriangleSetTopologyModifier, TopologicalChangeManager or TriangleCollisionModel not found";
			sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
			return;
		}
	}
}
