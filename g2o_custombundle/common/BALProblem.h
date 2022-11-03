#ifndef BALPROBLEM_H
#define BALPROBLEM_H

#include <stdio.h>
#include <string>
#include <iostream>

//传入数据文件（包含原始点，相机位姿，相机与点之间的关系）并解析。
//归一化，增加噪声的处理
class BALProblem
{
public:
    explicit BALProblem(const std::string& filename, bool use_quaternions = false);
    ~BALProblem(){
        delete[] point_index_;
        delete[] camera_index_;
        delete[] observations_;
        delete[] parameters_;
    }

    void WriteToFile(const std::string& filename)const;
    void WriteToPLYFile(const std::string& filename)const;

    void Normalize();

    void Perturb(const double rotation_sigma,
                 const double translation_sigma,
                 const double point_sigma);
    
    
    int camera_block_size()             const{ return use_quaternions_? 10 : 9;  }
    int point_block_size()              const{ return 3;                         }             
    int num_cameras()                   const{ return num_cameras_;              }
    int num_points()                    const{ return num_points_;               }
    int num_observations()              const{ return num_observations_;         }
    int num_parameters()                const{ return num_parameters_;           }
    const int* point_index()            const{ return point_index_;              }
    const int* camera_index()           const{ return camera_index_;             }
    const double* observations()        const{ return observations_;             }
    const double* parameters()          const{ return parameters_;               }
    const double* cameras()             const{ return parameters_;               }
    const double* points()              const{ return parameters_ + camera_block_size() * num_cameras_; }
    double* mutable_cameras()                { return parameters_;               }
    double* mutable_points()                 { return parameters_ + camera_block_size() * num_cameras_; }    //跳过了cameras数据，返回point数据的第一个

    double* mutable_camera_for_observation(int i){
        return mutable_cameras() + camera_index_[i] * camera_block_size();
    }

    double* mutable_point_for_observation(int i){
        return mutable_points() + point_index_[i] * point_block_size();
    }

    const double* camera_for_observation(int i)const {
        return cameras() + camera_index_[i] * camera_block_size();
    }

    const double* point_for_observation(int i)const {
        return points() + point_index_[i] * point_block_size();
    }


private:
    void CameraToAngelAxisAndCenter(const double* camera,
                                    double* angle_axis,
                                    double* center)const;

    void AngleAxisAndCenterToCamera(const double* angle_axis,
                                    const double* center,
                                    double* camera)const;

    int num_cameras_;
    int num_points_;
    int num_observations_;
    int num_parameters_;
    bool use_quaternions_;

    int* point_index_;
    int* camera_index_;
    double* observations_;
    double* parameters_;                             //inital.ply中最后面的一坨数据分为，上部分是9维的camera数据，旋转（轴角或四元数），平移，相机部分内参（f,K1,K2）

};

#endif // BALProblem.h
