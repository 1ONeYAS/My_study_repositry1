#include "BALProblem.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>


#include "tools/random.h"
#include "tools/rotation.h"


typedef Eigen::Map<Eigen::VectorXd> VectorRef;
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;

template<typename T>
void FscanfOrDie(FILE *fptr, const char *format, T *value){
        int num_scanned = fscanf(fptr, format, value);
        if(num_scanned != 1)
            std::cerr<< "Invalid UW data file. ";
}

void PerturbPoint3(const double sigma, double* point)
{
  for(int i = 0; i < 3; ++i)
    point[i] += RandNormal()*sigma;
}

//返回容器中从大到小排列出来在中间的数（容器中的中位数）
double Median(std::vector<double>* data){
  int n = data->size();
  std::vector<double>::iterator mid_point = data->begin() + n/2;
  std::nth_element(data->begin(),mid_point,data->end());
  return *mid_point;
}

//数据处理，将file中的文件数据填入相应的变量中（是否用寺院是表示由一个bool变量控制）
BALProblem::BALProblem(const std::string& filename, bool use_quaternions){
  FILE* fptr = fopen(filename.c_str(), "r");

    
  if (fptr == NULL) {
    std::cerr << "Error: unable to open file " << filename;
    return;
  };

  // This wil die horribly on invalid files. Them's the breaks.
  FscanfOrDie(fptr, "%d", &num_cameras_);              //总共有的相机个数（总共有多少帧）
  FscanfOrDie(fptr, "%d", &num_points_);                  //获取总共由多少关键点
  FscanfOrDie(fptr, "%d", &num_observations_);     //获取总共由多少个关系（总的观测数）

  std::cout << "Header: " << num_cameras_
            << " " << num_points_
            << " " << num_observations_;

  point_index_ = new int[num_observations_];                             //点索引
  camera_index_ = new int[num_observations_];                         //相机索引
  observations_ = new double[2 * num_observations_];

  num_parameters_ = 9 * num_cameras_ + 3 * num_points_;   //9维位姿，3维坐标点，总共需要多大的数组
  parameters_ = new double[num_parameters_];

  for (int i = 0; i < num_observations_; ++i) {
    FscanfOrDie(fptr, "%d", camera_index_ + i);
    FscanfOrDie(fptr, "%d", point_index_ + i);
    for (int j = 0; j < 2; ++j) {
      FscanfOrDie(fptr, "%lf", observations_ + 2*i + j);
    }
  }
 //获取文本中剩下的需要获取的参数（主要就是,9维camera，以及剩下的3维点坐标。）
  for (int i = 0; i < num_parameters_; ++i) {
    FscanfOrDie(fptr, "%lf", parameters_ + i);
  }

  fclose(fptr);
  //原来的数据结构，轴角【3】-平移【3】-f、k1、k2【3】
//当用四元数表示旋转的时候，从轴角的三维，变到四维。参数由9变成了10
//下面用两个数据指针来给9维的轴角坐标在第3个后面后面加上一维到10维【变成4】
  use_quaternions_ = use_quaternions;
  if (use_quaternions) {
    // Switch the angle-axis rotations to quaternions.
    num_parameters_ = 10 * num_cameras_ + 3 * num_points_;
    double* quaternion_parameters = new double[num_parameters_];
    double* original_cursor = parameters_;
    double* quaternion_cursor = quaternion_parameters;
    for (int i = 0; i < num_cameras_; ++i) {
      AngleAxisToQuaternion(original_cursor, quaternion_cursor);
      quaternion_cursor += 4;
      original_cursor += 3;
      for (int j = 4; j < 10; ++j) {
       *quaternion_cursor++ = *original_cursor++;
      }
    }
    // Copy the rest of the points.
    for (int i = 0; i < 3 * num_points_; ++i) {
      *quaternion_cursor++ = *original_cursor++;
    }
    // Swap in the quaternion parameters.
    delete []parameters_;
    parameters_ = quaternion_parameters;
  }
}

//没有研究
void BALProblem::WriteToFile(const std::string& filename)const{
  FILE* fptr = fopen(filename.c_str(),"w");
  
  if(fptr == NULL)
  {
    std::cerr<<"Error: unable to open file "<< filename;
    return;
  }

  fprintf(fptr, "%d %d %d %d\n", num_cameras_, num_cameras_, num_points_, num_observations_);

  for(int i = 0; i < num_observations_; ++i){
    fprintf(fptr, "%d %d", camera_index_[i], point_index_[i]);
    for(int j = 0; j < 2; ++j){
      fprintf(fptr, " %g", observations_[2*i + j]);
    }
    fprintf(fptr,"\n");
  }

  for(int i = 0; i < num_cameras(); ++i)
  {
    double angleaxis[9];
    if(use_quaternions_){
      //OutPut in angle-axis format.
      QuaternionToAngleAxis(parameters_ + 10 * i, angleaxis);
      memcpy(angleaxis + 3, parameters_ + 10 * i + 4, 6 * sizeof(double));
    }else{
      memcpy(angleaxis, parameters_ + 9 * i, 9 * sizeof(double));
    }
    for(int j = 0; j < 9; ++j)
    {
      fprintf(fptr, "%.16g\n",angleaxis[j]);
    }
  }

  const double* points = parameters_ + camera_block_size() * num_cameras_;
  for(int i = 0; i < num_points(); ++i){
    const double* point = points + i * point_block_size();
    for(int j = 0; j < point_block_size(); ++j){
      fprintf(fptr,"%.16g\n",point[j]);
    }
  }

  fclose(fptr);
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::WriteToPLYFile(const std::string& filename)const{
  std::ofstream of(filename.c_str());

  of<< "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << num_cameras_ + num_points_
    << '\n' << "property float x"
    << '\n' << "property float y"
    << '\n' << "property float z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;

    //从9维数据中提取除平移【3】的三维（即相机坐标离原点的距离）（并在ply文件中显示出来）
    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for(int i = 0; i < num_cameras(); ++i){
      const double* camera = cameras() + camera_block_size() * i;    //获获取parameters_（相机与点关系后面的数据）中上边，有关相机9维的数据（下面是点的3维坐标信息）
      CameraToAngelAxisAndCenter(camera, angle_axis, center);
      of << center[0] << ' ' << center[1] << ' ' << center[2]<<' ' 
         << "0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double* points = parameters_ + camera_block_size() * num_cameras_;
    for(int i = 0; i < num_points(); ++i){                                                             //对于初始ply文件中剩下的所有点（坐标点）进行处理
      const double* point = points + i * point_block_size();
      for(int j = 0; j < point_block_size(); ++j){                                               //第二个循环按顺序，将point位姿点进行处理
        of << point[j] << ' ';
      }
      of << "255 255 255\n";
    }
    of.close();
}

//将9维数据进行分割，分割除轴角（或者四元数），平移向量（相机中心）
void BALProblem::CameraToAngelAxisAndCenter(const double* camera, 
                                            double* angle_axis,
                                            double* center) const{
    VectorRef angle_axis_ref(angle_axis,3);    //将没填充的向量进行'内存'映射绑定
    if(use_quaternions_){
      QuaternionToAngleAxis(camera, angle_axis);    //如果是四元数，就通过四元数转换轴角函数，直接求出轴角
    }else{
      angle_axis_ref = ConstVectorRef(camera,3);  //如果是轴角，就通过简单的映射9维的前3维。
    }

    // c = -R't
    Eigen::VectorXd inverse_rotation = -angle_axis_ref;
    AngleAxisRotatePoint(inverse_rotation.data(),     //对平移向量（相机中心）【这里像是对平移进行了旋转】？？？？？---希望在平移图中知道是怎么计算的
                         camera + camera_block_size() - 6,
                         center);
    VectorRef(center,3) *= -1.0;
}

//主要用在加噪声或者正态化后，将轴角，平移数据放回9维的数据中
void BALProblem::AngleAxisAndCenterToCamera(const double* angle_axis,
                                            const double* center,
                                            double* camera) const{
    ConstVectorRef angle_axis_ref(angle_axis,3);
    if(use_quaternions_){
      AngleAxisToQuaternion(angle_axis,camera);
    }
    else{
      VectorRef(camera, 3) = angle_axis_ref;
    }

    // t = -R * c 
    AngleAxisRotatePoint(angle_axis,center,camera+camera_block_size() - 6);
    VectorRef(camera + camera_block_size() - 6,3) *= -1.0;
}


//数据标准化
void BALProblem::Normalize(){
  // Compute the marginal median of the geometry
  std::vector<double> tmp(num_points_);      //指定元素个数进行vector的初始化
  Eigen::Vector3d median;
  double* points = mutable_points();                   //temp按顺序将所有位置第一个放在最上面
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < num_points_; ++j){
      tmp[j] = points[3 * j + i];      
    }
    median(i) = Median(&tmp);    //分别将所有点的X，Y，Z坐标中的中位数找出来，（X的中位数为第一个）
  }

  for(int i = 0; i < num_points_; ++i){
    VectorRef point(points + 3 * i, 3);      //？每次映射3个点为一个向量vector3d的类型
    tmp[i] = (point - median).lpNorm<1>();     //？每个关键点与中位数点的3个坐标的差并取1范数
  }

  const double median_absolute_deviation = Median(&tmp);   //所有关键点与中位数差1范数的“中位数”

  // Scale so that the median absolute deviation of the resulting
  // reconstruction is 100

  const double scale = 100.0 / median_absolute_deviation;

  // X = scale * (X - median)
  for(int i = 0; i < num_points_; ++i){
    VectorRef point(points + 3 * i, 3);
    point = scale * (point - median);
  }

//位姿参数中，平移减去中位数并乘以比例系数，并转化回相机参数类型
  double* cameras = mutable_cameras();
  double angle_axis[3];
  double center[3];
  for(int i = 0; i < num_cameras_ ; ++i){
    double* camera = cameras + camera_block_size() * i;
    CameraToAngelAxisAndCenter(camera, angle_axis, center);
    // center = scale * (center - median)
    VectorRef(center,3) = scale * (VectorRef(center,3)-median);
    AngleAxisAndCenterToCamera(angle_axis, center,camera);
  }
}


//对旋转，平移，坐标点，增加噪声
void BALProblem::Perturb(const double rotation_sigma, 
                         const double translation_sigma,
                         const double point_sigma){
   assert(point_sigma >= 0.0);                  //assert判断是否满足条件，如果不满足，代码终止并返回错误内容（像单片机标准库中的assert一样）
   assert(rotation_sigma >= 0.0);
   assert(translation_sigma >= 0.0);

   double* points = mutable_points();
   if(point_sigma > 0){
     for(int i = 0; i < num_points_; ++i){      //对每个点增加扰动（随机数乘以方差的扰动）
       PerturbPoint3(point_sigma, points + 3 * i);
     }
   }

   for(int i = 0; i < num_cameras_; ++i){
     double* camera = mutable_cameras() + camera_block_size() * i;    //找到循环中每个相机数据的第一个

     double angle_axis[3];
     double center[3];
     // Perturb in the rotation of the camera in the angle-axis
     // representation
     CameraToAngelAxisAndCenter(camera, angle_axis, center);
     if(rotation_sigma > 0.0){
       PerturbPoint3(rotation_sigma, angle_axis);       
     }
     AngleAxisAndCenterToCamera(angle_axis, center,camera);    //对轴角增加扰动，并返回camera*

     if(translation_sigma > 0.0)
        PerturbPoint3(translation_sigma, camera + camera_block_size() - 6);  //对平移增加扰动直接在camera*中进行
   }
}