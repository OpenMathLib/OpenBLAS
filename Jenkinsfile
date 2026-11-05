pipeline {
    agent { 
        docker {
            image 'osuosl/ubuntu-s390x'
        }
    }
    stages {
        stage('Build') {
            steps {
                sh 'make MAKE_NB_JOBS=2 clean && make MAKE_NB_JOBS=2'
            }
        }
        stage('CMakeBuild') {
            steps {
                sh 'sudo apt update && sudo apt install cmake -y && make -j 2 clean && rm -rf build && mkdir build && cd build && cmake -DDYNAMIC_ARCH=1 .. && make -j 2'
            }
        }
    }
}
