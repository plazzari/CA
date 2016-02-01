#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/features/normal_3d.h>
#include <pcl/surface/gp3.h>
#include <pcl/io/vtk_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/sac_model_normal_plane.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/extract_clusters.h>
#include <iostream>
#include <fstream>

int
main (int argc, char** argv)
{
  pcl::PointCloud<pcl::PointXYZ>:: Ptr cloud_p (new pcl::PointCloud<pcl::PointXYZ>), cloud_f (new pcl::PointCloud<pcl::PointXYZ>);
 pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal> seg;

  // Load input file into a PointCloud<T> with an appropriate type
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>), cloudVP (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PCLPointCloud2 cloud_blob, cloud_blobVP;
  pcl::io::loadPCDFile (argv[1], cloud_blob);
  pcl::fromPCLPointCloud2 (cloud_blob, *cloud);

  pcl::io::loadPCDFile (argv[2], cloud_blobVP); // Load viewpoint dataset
  pcl::fromPCLPointCloud2 (cloud_blobVP, *cloudVP);
  //* the data should be available in cloud

  // Normal estimation*
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>), normals_p (new pcl::PointCloud<pcl::Normal>), normals_f (new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud (cloud);
  n.setSearchMethod (tree);
  n.setKSearch (20);
  n.compute (*normals);
  //* normals should not contain the point normals + surface curvatures

  // Create the filtering object to flip normals
  pcl::ExtractIndices<pcl::PointXYZ> ex_to_flip;
  pcl::ExtractIndices<pcl::Normal> ex_n_to_flip;

  int j = 0, nr_normals = (int) normals->points.size ();
  std::cerr << "Number of normals " << nr_normals << std::endl;
  while (j < nr_normals) 
  {

  float p_x = (float) cloud->points[j].x; float vp_x = (float) cloudVP->points[j].x;
  float p_y = (float) cloud->points[j].y; float vp_y = (float) cloudVP->points[j].y;
  float p_z = (float) cloud->points[j].z; float vp_z = (float) cloudVP->points[j].z;

  float n_x = (float) normals->points[j].normal[0]; 
  float n_y = (float) normals->points[j].normal[1]; 
  float n_z = (float) normals->points[j].normal[2]; 
  float v_v = n_x * ( vp_x -p_x) + n_y * ( vp_y -p_y) + n_z * ( vp_z -p_z);
  if (v_v > 0)
  {
  normals->points[j].normal[0] = - normals->points[j].normal[0];
  normals->points[j].normal[1] = - normals->points[j].normal[1];
  normals->points[j].normal[2] = - normals->points[j].normal[2];
  }
  j++;
  }

  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals_p (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);
  //* cloud_with_normals = cloud + normals

  // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
  pcl::PolygonMesh triangles;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (10); //  0.025

  // Set typical values for the parameters
  gp3.setMu (2.5);
  gp3.setMaximumNearestNeighbors (100);
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (triangles);

  // Additional vertex information
  std::vector<int> parts = gp3.getPartIDs();
  std::vector<int> states = gp3.getPointStates();

// Create the segmentation object for the planar model and set all the parameters
  seg.setOptimizeCoefficients (true);
  seg.setModelType (pcl::SACMODEL_NORMAL_PLANE);
  seg.setNormalDistanceWeight (0.1);
  seg.setMethodType (pcl::SAC_RANSAC);
  seg.setMaxIterations (1000);
  seg.setDistanceThreshold (0.1); // 0.1


  pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients ());
  pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());

  // Create the filtering object
  pcl::ExtractIndices<pcl::PointXYZ> extract;
  pcl::ExtractIndices<pcl::Normal> extract_n;

  // Creating the KdTree object for the search method of the cluster extraction

   pcl::search::KdTree<pcl::PointXYZ>::Ptr tree3 (new pcl::search::KdTree<pcl::PointXYZ>);
   pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZ>);
   pcl::PointCloud<pcl::Normal>::Ptr normals_filtered (new pcl::PointCloud<pcl::Normal>);

///////////////////////////////////////////////////////// 
// START DOING THE JOB !!!
///////////////////////////////////////////////////////// 
  int ij =0;
  int i = 0, nr_points = (int) cloud->points.size ();
  // While 30% of the original cloud is still there
  while (cloud->points.size () > 0.02 * nr_points) //0.3
  {
    // Segment the largest planar component from the remaining cloud
    seg.setInputCloud (cloud);
    seg.setInputNormals (normals);
    seg.segment (*inliers, *coefficients);
    if (inliers->indices.size () == 0)
    {
      std::cerr << "Could not estimate a planar model for the given dataset." << std::endl;
      break;
    }

    // Extract the inliers
    extract.setInputCloud (cloud);
    extract.setIndices (inliers);
    extract.setNegative (false);
    extract.filter (*cloud_p);

    extract_n.setInputCloud (normals);
    extract_n.setIndices (inliers);
    extract_n.setNegative (false);
    extract_n.filter (*normals_p);
    std::cerr << "PointCloud representing the planar component: " << cloud_p->width * cloud_p->height << " data points." << std::endl;


    // Create the filtering object
    extract.setNegative (true);
    extract.filter (*cloud_f);
    cloud.swap (cloud_f);

    // Extract the inliers for normals
    extract_n.setNegative (true);
    extract_n.filter (*normals_f);
    normals.swap (normals_f);
   
    *cloud_filtered=*cloud_p;
    *normals_filtered=*normals_p;
    tree3->setInputCloud (cloud_filtered);

    std::vector<pcl::PointIndices> cluster_indices;
    pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
    ec.setClusterTolerance (3); // 2cm
    ec.setMinClusterSize (5);
    ec.setMaxClusterSize (25000);
    ec.setSearchMethod (tree3);
    ec.setInputCloud (cloud_filtered);
    ec.extract (cluster_indices);
    int j = 0;
    for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin (); it != cluster_indices.end (); ++it)
    {
      pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster (new pcl::PointCloud<pcl::PointXYZ>);
      pcl::PointCloud<pcl::Normal>::Ptr normals_cluster (new pcl::PointCloud<pcl::Normal>);
      pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals_cluster (new pcl::PointCloud<pcl::PointNormal>);

      for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
      cloud_cluster->points.push_back (cloud_filtered->points[*pit]); //*
      for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit)
      normals_cluster->points.push_back (normals_filtered->points[*pit]); //*
      pcl::concatenateFields (*cloud_cluster, *normals_cluster, *cloud_with_normals_cluster);
      cloud_with_normals_cluster->width = cloud_with_normals_cluster->points.size ();
      cloud_with_normals_cluster->height = 1;
      cloud_with_normals_cluster->is_dense = true;

      std::cout << "PointCloud representing the Cluster:" << ij << " " << cloud_cluster->points.size () << " data points." << std::endl;


// write pcd file for the specific face
      std::stringstream ss;
      ss << argv[3] << ij << ".pcd";
      pcl::io::savePCDFile(ss.str (), *cloud_with_normals_cluster, false);	

// write ascii file for the plane face

      std::stringstream ss1;
      std::ofstream myfile;
      ss1 << argv[3] << ij << ".coeff.txt";
      std:: string fileName = ss1.str();

      myfile.open(fileName.c_str());

      if (myfile.is_open())
      {
       myfile << coefficients->values[0] << " ";
       myfile << coefficients->values[1] << " ";
       myfile << coefficients->values[2] << " "; 
       myfile << coefficients->values[3] << " ";
       myfile.close();
      }

      j++;
      ij++;
    }

    i++;
  }
  // Finish
  return (0);
}
