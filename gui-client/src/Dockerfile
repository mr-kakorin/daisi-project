# Сборка ---------------------------------------

# В качестве базового образа для сборки используем gcc:latest
FROM gcc:latest as build

# Установим рабочую директорию для сборки GoogleTest
#WORKDIR /build

# Скачаем все необходимые пакеты и выполним сборку GoogleTest
# Такая длинная команда обусловлена тем, что
# Docker на каждый RUN порождает отдельный слой,
# Влекущий за собой, в данном случае, ненужный оверхед
RUN apt-get update && \
    apt-get install -y \
      libboost-dev libboost-program-options-dev \
      libgtest-dev \
      cmake \
      gconf-service libasound2 libatk1.0-0 libc6 libcairo2 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libgcc1 libgconf-2-4 libgdk-pixbuf2.0-0 libglib2.0-0 libgtk-3-0 libnspr4 libpango-1.0-0 libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 libxcomposite1 libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 libxss1 libxtst6 ca-certificates fonts-liberation libappindicator1 libnss3 lsb-release xdg-utils wget \
    && \
    #cmake -DCMAKE_BUILD_TYPE=Debug /usr/src/gtest && \
    #cmake --build . && \
    #mv lib*.a /usr/lib && \
    apt-get -y install libfontconfig1 && \
    apt-get -y install mesa-common-dev \
    libglu1-mesa-dev -y  \
    libxcb1-dev libicu-dev libssl-dev \
    #linux-headers-generic curl git libxrender-dev \
    #libpng3 libpng-dev libjpeg-turbo8 libjpeg-turbo8-dev libicu55 libicu-dev \
    libgles2-mesa libgles2-mesa-dev && \
    #libfreetype6-dev libsqlite3-dev \
    #libgstreamer0.10-0 libgstreamer0.10-dev \
    #libogg-dev libvorbis-dev bzip2 gperf bison ruby flex && \
    curl -sSL https://download.qt.io/official_releases/qt/5.9/5.9.0/submodules/qtbase-opensource-src-5.9.0.tar.xz | tar xJ && \
    #ls && \
    #cd qtbase-opensource-src-5.9.0 && \
    #bash ./configure --help && \
    #bash ./configure -opensource -confirm-license -static -no-accessibility -no-qml-debug \
    #   -no-harfbuzz -qt-pcre -no-dbus -nomake tools -no-openssl\
    #   -no-xkbcommon-evdev -no-xinput2 -no-xcb-xlib -no-glib -qt-xcb -no-compile-examples -nomake examples \
    #   -no-gif -qt-doubleconversion -no-gtk -platform linux-g++-64 && \
    #make install && \
    apt-get -y install gfortran gfortran-6 && \
    update-alternatives --install /usr/bin/gfortran gfortran /usr/local/bin/gfortran 999 && \
    update-alternatives --install /usr/bin/gfortran-6 gfortran /usr/local/bin/gfortran-6 999 && \
    apt-get -y install libblas-dev && \
    apt-get -y install libopenblas-dev \
    apt-get -y install liblapack-dev && \
    apt-get -y install libboost-dev && \
    apt-get -y install libarmadillo-dev && \
    apt-get -y install libvtk6.2

ENV PATH $QT_DIST/bin:$PATH

# Скопируем директорию /src в контейнер
ADD ./ /app/src

# Установим рабочую директорию для сборки проекта
WORKDIR /app/build

# Выполним сборку нашего проекта, а также его тестирование
RUN cmake ../src && \
    cmake --build .

# Запуск ---------------------------------------

# В качестве базового образа используем ubuntu:latest
FROM ubuntu:latest

# Добавим пользователя, потому как в Docker по умолчанию используется root
# Запускать незнакомое приложение под root'ом неприлично :)
RUN groupadd -r sample && useradd -r -g sample sample
USER sample

# Установим рабочую директорию нашего приложения
WORKDIR /app

# Скопируем приложение со сборочного контейнера в рабочую директорию
COPY --from=build /app/build/daisi_client .

# Установим точку входа
ENTRYPOINT ["./daisi_client"]