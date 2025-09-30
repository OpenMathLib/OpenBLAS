pipeline {
    agent { 
        docker {
            image 'osuosl/ubuntu-s390x'
        }
    }
    stages {
        stage('Build') {
            steps {
                sh 'make clean && make'
            }
        }
        stage('CMakeBuild') {
            steps {
                sh 'sudo apt update && sudo apt install cmake -y && make clean && rm -rf build && mkdir build && cd build && cmake -DDYNAMIC_ARCH=1 .. && make'
            }
        }
    }
}
