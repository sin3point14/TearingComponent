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

		Size nbTriangles = this->m_topology->getNbTriangles();

    	const auto& triangleInf = this->triangleInfo.getValue();

		sofa::component::forcefield::TriangularFEMForceField<DataTypes>::addForce(mparams, f, x, v);
		const VecCoord& x1 = this->mstate->read(core::ConstVecCoordId::position())->getValue();
        std::vector<sofa::type::Vector3> vertices;
        for(Size i=0; i< nbTriangles; ++i)
        {
            const sofa::core::topology::BaseMeshTopology::Triangle& tri = this->m_topology->getTriangle(i);
            Index a = tri[0];
            Index b = tri[1];
            Index c = tri[2];
            Coord center = (x1[a]+x1[b]+x1[c])/3;
            Coord d = triangleInf[i].principalStressDirection*2.5; //was 0.25
            vertices.push_back(sofa::type::Vector3(center));
            vertices.push_back(sofa::type::Vector3(center+d));
        }
        // vparams->drawTool()->drawLines(vertices, 1, sofa::type::RGBAColor(1, 0, 1, 1));
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

		if (m_triangleMod == nullptr || m_triangleGeo == nullptr)
		{
			msg_error() << "No TriangleSetTopologyModifier or TriangleSetGeometryAlgorithms component found";
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
