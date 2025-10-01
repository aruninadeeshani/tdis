// Copyright (C) 2025 Dmitry Romanov

#pragma once

#include <JANA/Components/JComponent.h>
#include <JANA/JApplication.h>
#include <JANA/Services/JServiceLocator.h>
#include <TGeoManager.h>
#include <spdlog/logger.h>

#include <memory>
#include <mutex>

#include "geometry/MtpcDetectorElement.h"
#include "logging/LogService.hpp"

namespace tdis::tracking {
    class TGeoGeometryService : public JService {
      public:
        explicit TGeoGeometryService() : JService() {}
        ~TGeoGeometryService() override = default;

        void Init() override;

        TGeoManager* GetGeoManager() const { return m_tgeo_manager; }

        TGeoVolume* GetTopVolume() const { return m_tgeo_manager->GetTopVolume(); }

        TGeoNode* GetTopNode() const { return m_tgeo_manager->GetTopNode(); }

        const std::vector<double>& GetPlanePositions() const { return m_plane_positions; }


      private:
        Parameter<std::string> m_tgeo_file{this, "acts:geometry", "g4sbs_mtpc.root","TGeo filename with geometry for ACTS"};
        Parameter<std::string> m_material_map_file{this, "acts:material_map", "", "JSON/CBOR material map file path"};
        Parameter<std::string> m_obj_output_file{this, "acts:output_obj", "", "Output file name to dump ACTS converted geometry as OBJ"};
        Parameter<std::string> m_ply_output_file{this, "acts:output_ply", "", "Output file name to dump ACTS converted geometry as PLY"};

        Service<tdis::services::LogService> m_service_log{this};

        // General acts log
        std::shared_ptr<spdlog::logger> m_log;

        /// Root TGeo Manater for TGeo Geometry
        TGeoManager* m_tgeo_manager = nullptr;

        // Plane positions
        std::vector<double> m_plane_positions;
    };

}   // namespace tdis::tracking