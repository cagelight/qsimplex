#include <meadow/noise.hh>

#include <random>
#include <source_location>

#include <QCoreApplication>
#include <QCommandLineParser>
#include <QImageWriter>
#include <QFile>
#include <QDebug>
#include <QMap>

std::mt19937_64 rng { std::random_device {} () };
std::uniform_real_distribution<double> rdist {-1024, 1024};
#define RNG rdist(rng)

struct RenderSettings {
	int dims;
	double scale = 32.0;
	int frames = 30;
	double anim = 0.1;
	std::function<void(QImage, QString suffix)> write_cb;
};

using RendFunc = std::function<void(RenderSettings const &)>;

// ================
// 2d
// ================

void rend_2d(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	for (size_t y = 0; y < rs.dims; y++) {
		double yv = (static_cast<double>(y) / rs.dims) * rs.scale;
		uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
		for (size_t x = 0; x < rs.dims; x++) {
			double xv = (static_cast<double>(x) / rs.dims) * rs.scale;
			line[x] = (noise.generate(xv, yv) + 1) * std::numeric_limits<uint16_t>::max() / 2;
		}
	}
	
	rs.write_cb(img, "");
}

// ================
// 3d
// ================


void rend_3d(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	double zr = RNG;
	
	for (size_t y = 0; y < rs.dims; y++) {
		double yv = (static_cast<double>(y) / rs.dims) * rs.scale;
		uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
		for (size_t x = 0; x < rs.dims; x++) {
			double xv = (static_cast<double>(x) / rs.dims) * rs.scale;
			line[x] = (noise.generate(xv, yv, zr) + 1) * std::numeric_limits<uint16_t>::max() / 2;
		}
	}
	
	rs.write_cb(img, "");
}

// ================
// 3da
// ================

void rend_3da(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	for (size_t i = 0; i < rs.frames; i++) {
		
		double zv = i * rs.anim;
		
		for (size_t y = 0; y < rs.dims; y++) {
			double yv = (static_cast<double>(y) / rs.dims) * rs.scale;
			uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
			for (size_t x = 0; x < rs.dims; x++) {
				double xv = (static_cast<double>(x) / rs.dims) * rs.scale;
				line[x] = (noise.generate(xv, yv, zv) + 1) * std::numeric_limits<uint16_t>::max() / 2;
			}
		}
		
		rs.write_cb(img, QStringLiteral("%1").arg(i, 5, 10, QLatin1Char('0')));
	}
}

// ================
// 4d
// ================

void rend_4d(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	double zr = RNG;
	double wr = RNG;
	
	for (size_t y = 0; y < rs.dims; y++) {
		double yv = (static_cast<double>(y) / rs.dims) * rs.scale;
		uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
		for (size_t x = 0; x < rs.dims; x++) {
			double xv = (static_cast<double>(x) / rs.dims) * rs.scale;
			line[x] = (noise.generate(xv, yv, zr, wr) + 1) * std::numeric_limits<uint16_t>::max() / 2;
		}
	}
	
	rs.write_cb(img, "");
}

// ================
// 4da
// ================

void rend_4da(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	for (size_t i = 0; i < rs.frames; i++) {
		
		double zv, wv;
		zv = wv = (i / static_cast<double>(rs.frames)) * meadow::brassica::pi_m<double>(2);
		zv = std::sin(zv) * rs.anim;
		wv = std::cos(wv) * rs.anim;
		
		for (size_t y = 0; y < rs.dims; y++) {
			double yv = (static_cast<double>(y) / rs.dims) * rs.scale;
			uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
			for (size_t x = 0; x < rs.dims; x++) {
				double xv = (static_cast<double>(x) / rs.dims) * rs.scale;
				line[x] = (noise.generate(xv, yv, zv, wv) + 1) * std::numeric_limits<uint16_t>::max() / 2;
			}
		}
		
		rs.write_cb(img, QStringLiteral("%1").arg(i, 5, 10, QLatin1Char('0')));
	}
}

// ================
// 4dt
// ================

void rend_4dt(RenderSettings const & rs) {
	QImage img { rs.dims, rs.dims, QImage::Format_Grayscale16 };
	meadow::noise::simplex noise;
	noise.seed();
	
	double TILE_M = rs.scale / meadow::brassica::pi_m<double>(2);
	
	for (size_t y = 0; y < rs.dims; y++) {
		double yv = (static_cast<double>(y) / rs.dims) * meadow::brassica::pi_m<double>(2);
		uint16_t * line = reinterpret_cast<uint16_t *>(img.scanLine(y));
		for (size_t x = 0; x < rs.dims; x++) {
			double xv = (static_cast<double>(x) / rs.dims) * meadow::brassica::pi_m<double>(2);
			line[x] = (noise.generate(std::sin(xv) * TILE_M, std::cos(xv) * TILE_M, std::sin(yv) * TILE_M, std::cos(yv) * TILE_M) + 1) * std::numeric_limits<uint16_t>::max() / 2;
		}
	}
	
	rs.write_cb(img, "");
}

QMap<QString, QPair<QString, RendFunc>> rend_funcs {
	{ QString("2d"),   { QString("2D algorithm"),                rend_2d  } },
	{ QString("3d"),   { QString("3D algorithm [default]"),      rend_3d  } },
	{ QString("3da"),  { QString("3D algorithm, non-loop anim"), rend_3da } },
	{ QString("4d"),   { QString("4D algorithm"),                rend_4d  } },
	{ QString("4da"),  { QString("4D algorithm, looping anim"),  rend_4da } },
	{ QString("4dt"),  { QString("4D algorithm, tiling"),        rend_4dt } }
};

QString generate_mode_help() {
	QString str;
	str.append("Render mode:\n");
	for (auto const & [k, v] : rend_funcs.toStdMap()) {
		str += QString("    %0: %1\n") .arg(k) .arg(v.first);
	}
	return str;
}

int main(int argc, char * * argv) {
	
	QCoreApplication app { argc, argv };
	
	QCommandLineParser parser;
	parser.addHelpOption();
	
	QCommandLineOption dimsOpt {QStringList() << "d" << "dims", "Dimensions (width and height) of the output image.", "dimensions"};
	parser.addOption(dimsOpt);
	QCommandLineOption scaleOpt {QStringList() << "s" << "scale", "Noise scale. (Default: 32)", "multiplier"};
	parser.addOption(scaleOpt);
	QCommandLineOption framesOpt {QStringList() << "f" << "frames", "Frame count, used only in animated modes. (Default: 30)", "count"};
	parser.addOption(framesOpt);
	QCommandLineOption animOpt {QStringList() << "a" << "anim", "Animation \"magnitude\", a catch-all to adjust animations. (Default: 0.1)", "mag"};
	parser.addOption(animOpt);
	QCommandLineOption modeOpt {QStringList() << "m" << "mode", generate_mode_help(), "mode"};
	parser.addOption(modeOpt);
	
	parser.addPositionalArgument("<output>", "Output path without extension.");
	
	parser.process(app);
	 
	if (!parser.positionalArguments().size()) {
		qDebug() << "output path is required!\n";
		parser.showHelp();
		exit(-1);
	}
	
	QString mode = "3d";
	RenderSettings rs;
	bool ok;
	
	// ================
	// DIMS
	// ================
	
	if (!parser.isSet(dimsOpt)) {
		qDebug() << "dimensions are required!\n";
		parser.showHelp();
		exit(-1);
	}
	
	rs.dims = parser.value(dimsOpt).toInt(&ok);
	if (!ok) {
		qDebug() << "invalid dimensions!\n";
		parser.showHelp();
		exit(-1);
	}
	
	if (rs.dims < 2) rs.dims = 2;
	if (rs.dims > 32768) rs.dims = 32768;
	
	// ================
	// SCALE
	// ================
	
	if (parser.isSet(scaleOpt)) {
		rs.scale = parser.value(scaleOpt).toDouble(&ok);
		if (!ok) {
			qDebug() << "invalid scale!\n";
			parser.showHelp();
			exit(-1);
		}
	}
	
	// ================
	// MODE
	// ================
	
	if (parser.isSet(modeOpt)) {
		mode = parser.value(modeOpt);
	}

	if (!rend_funcs.contains(mode)) {
		qDebug() << "invalid mode!\n";
		parser.showHelp();
		exit(-1);
	}
	
	// ================
	// FRAMES
	// ================
	
	if (parser.isSet(framesOpt)) {
		rs.frames = parser.value(framesOpt).toInt(&ok);
		if (!ok) {
			qDebug() << "invalid frames value!\n";
			parser.showHelp();
			exit(-1);
		}
	}
	
	// ================
	// ANIM
	// ================
	
	if (parser.isSet(animOpt)) {
		rs.anim = parser.value(animOpt).toDouble(&ok);
		if (!ok) {
			qDebug() << "invalid anim value!\n";
			parser.showHelp();
			exit(-1);
		}
	}
	
	// ================
	
	rs.write_cb = [&](QImage img, QString suffix){
		QString out_path = parser.positionalArguments().at(0) + (suffix.isEmpty() ? ".png" : "_" + suffix + ".png");
		QFile out { out_path };
		if (!out.open(QIODevice::WriteOnly)) {
			qDebug() << "output could not be opened for writing!\n";
			exit(-1);
		}
		
		QImageWriter imgw { &out, "png" };
		imgw.write(img);
	};
	
	rend_funcs[mode].second(rs);
	
	return 0;
}
