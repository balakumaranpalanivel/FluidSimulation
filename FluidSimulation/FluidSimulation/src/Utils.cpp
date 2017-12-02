#include "Utils.h"

namespace Utils
{
	void DrawBillboard(CCamera *camera, GLuint *tex, glm::vec3 p, float width)
	{
		float hw = 0.5*width;
		glm::vec3 look = glm::normalize(camera->position - p);
		glm::vec3 right = glm::normalize(glm::cross(camera->up, look));
		glm::vec3 up = glm::cross(look, right);

		glm::mat4 mat = glm::transpose(glm::mat4(right.x, up.x, look.x, p.x,
			right.y, up.y, look.y, p.y,
			right.z, up.z, look.z, p.z,
			0.0, 0.0, 0.0, 1.0));
		glPushMatrix();
		glMultMatrixf((GLfloat*)&mat);

		glEnable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glBindTexture(GL_TEXTURE_2D, tex[0]);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(-hw, -hw, 0.0f);
		glTexCoord2f(1.0f, 0.0f); glVertex3f(hw, -hw, 0.0f);
		glTexCoord2f(1.0f, 1.0f); glVertex3f(hw, hw, 0.0f);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(-hw, hw, 0.0f);
		glEnd();
		glEnable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glPopMatrix();
	}

	void DrawGrid()
	{
		// draw axis'
		float len = 10.0;
		glLineWidth(3.0);
		glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);   // x
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(len, 0.0, 0.0);
		glColor3f(0.0, 1.0, 0.0);   // y
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, len, 0.0);
		glColor3f(0.0, 0.0, 1.0);   // z
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, len);
		glEnd();

		// draw outline around xy, zy planes
		glLineWidth(2.0);
		glColor4f(0.0, 0.0, 0.0, 0.3);
		glBegin(GL_LINES);
		glVertex3f(0.0, len, 0.0);
		glVertex3f(len, len, 0.0);
		glVertex3f(len, len, 0.0);
		glVertex3f(len, 0.0, 0.0);
		glVertex3f(0.0, len, 0.0);
		glVertex3f(0.0, len, len);
		glVertex3f(0.0, len, len);
		glVertex3f(0.0, 0.0, len);
		glEnd();


		// draw xz plane grid
		float spacing = 0.25;
		int yLines = 120;
		int zLines = 120;
		float height = (float)yLines * spacing;
		float width = (float)zLines * spacing;

		float z = spacing;
		glLineWidth(1.0);
		glColor4f(0.0, 0.0, 0.0, 0.2);
		glBegin(GL_LINES);
		for (int i = 0; i < yLines; i++) {
			glVertex3f(0.0, 0.0, z);
			glVertex3f(width, 0.0, z);
			z += spacing;
		}

		float x = spacing;
		for (int i = 0; i < zLines; i++) {
			glVertex3f(x, 0.0, 0.0);
			glVertex3f(x, 0.0, height);
			x += spacing;
		}
		glEnd();

	}

	void DrawWireFrameCube(glm::vec3 pos, float size)
	{
		float h = 0.5*size;
		glBegin(GL_LINES);
		glVertex3f(pos.x - h, pos.y - h, pos.z - h);
		glVertex3f(pos.x + h, pos.y - h, pos.z - h);
		glVertex3f(pos.x - h, pos.y - h, pos.z - h);
		glVertex3f(pos.x - h, pos.y + h, pos.z - h);
		glVertex3f(pos.x - h, pos.y - h, pos.z - h);
		glVertex3f(pos.x - h, pos.y - h, pos.z + h);

		glVertex3f(pos.x + h, pos.y + h, pos.z + h);
		glVertex3f(pos.x - h, pos.y + h, pos.z + h);
		glVertex3f(pos.x + h, pos.y + h, pos.z + h);
		glVertex3f(pos.x + h, pos.y - h, pos.z + h);
		glVertex3f(pos.x + h, pos.y + h, pos.z + h);
		glVertex3f(pos.x + h, pos.y + h, pos.z - h);

		glVertex3f(pos.x - h, pos.y + h, pos.z + h);
		glVertex3f(pos.x - h, pos.y - h, pos.z + h);
		glVertex3f(pos.x - h, pos.y + h, pos.z + h);
		glVertex3f(pos.x - h, pos.y + h, pos.z - h);

		glVertex3f(pos.x + h, pos.y - h, pos.z + h);
		glVertex3f(pos.x - h, pos.y - h, pos.z + h);
		glVertex3f(pos.x + h, pos.y - h, pos.z + h);
		glVertex3f(pos.x + h, pos.y - h, pos.z - h);

		glVertex3f(pos.x + h, pos.y + h, pos.z - h);
		glVertex3f(pos.x + h, pos.y - h, pos.z - h);
		glVertex3f(pos.x + h, pos.y + h, pos.z - h);
		glVertex3f(pos.x - h, pos.y + h, pos.z - h);

		glEnd();

	}

	void DrawWireFrameCube(glm::vec3 pos, float width, float height, float depth)
	{
		float hw = 0.5*width;
		float hh = 0.5*height;
		float hd = 0.5*depth;

		glBegin(GL_LINES);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z - hd);
		glVertex3f(pos.x + hw, pos.y - hh, pos.z - hd);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z - hd);
		glVertex3f(pos.x - hw, pos.y + hh, pos.z - hd);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z - hd);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z + hd);

		glVertex3f(pos.x + hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x - hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y - hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y + hh, pos.z - hd);

		glVertex3f(pos.x - hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z + hd);
		glVertex3f(pos.x - hw, pos.y + hh, pos.z + hd);
		glVertex3f(pos.x - hw, pos.y + hh, pos.z - hd);

		glVertex3f(pos.x + hw, pos.y - hh, pos.z + hd);
		glVertex3f(pos.x - hw, pos.y - hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y - hh, pos.z + hd);
		glVertex3f(pos.x + hw, pos.y - hh, pos.z - hd);

		glVertex3f(pos.x + hw, pos.y + hh, pos.z - hd);
		glVertex3f(pos.x + hw, pos.y - hh, pos.z - hd);
		glVertex3f(pos.x + hw, pos.y + hh, pos.z - hd);
		glVertex3f(pos.x - hw, pos.y + hh, pos.z - hd);

		glEnd();
	}

	float Lerp(float x1, float x2, float t)
	{
		return x1 + t*(x2 - x1);
	}

	float SmoothStep(float t)
	{
		return t*t*(3 - 2 * t);
	}

	std::vector<glm::vec3> CreatePointPanel(float width, float height,
		float spacing, int numLayers,
		glm::vec3 w, glm::vec3 h, bool isStaggered)
	{
		// adjust spacing to fit width/height constraints
		int spaces = floor(width / spacing);
		float xpad = width / spaces;
		int nx = spaces + 1;

		spaces = floor(height / spacing);
		float ypad = height / spaces;
		int ny = spaces + 1;

		float zpad = spacing;
		int nz = numLayers;

		w = glm::normalize(w);
		h = glm::normalize(h);
		glm::vec3 normal = glm::cross(w, h);

		std::vector<glm::vec3> points;
		glm::vec3 start = (float)(-0.5*width)*w -
			(float)(0.5*height)*h -
			(float)(0.5*(nz - 1)*zpad)*normal;
		for (int k = 0; k<nz; k++) {
			for (int j = 0; j<ny; j++) {
				for (int i = 0; i<nx; i++) {
					glm::vec3 p = start + i*xpad*w + j*ypad*h + k*zpad*normal;
					points.push_back(p);

					if (isStaggered && i != nx - 1 && j != ny - 1 && (k != nz - 1 || nz == 1)) {
						p = p + (float)(0.5*xpad)*w +
							(float)(0.5*ypad)*h +
							(float)(0.5*zpad)*normal;
						points.push_back(p);
					}
				}
			}
		}

		return points;

	}

	std::vector<glm::vec3> TranslatePoints(std::vector<glm::vec3> points, glm::vec3 trans)
	{
		std::vector<glm::vec3> newPoints;
		for (unsigned int i = 0; i<points.size(); i++) {
			newPoints.push_back(points[i] + trans);
		}
		return newPoints;
	}

	std::vector<glm::vec3> MergePoints(std::vector<glm::vec3> points1,
		std::vector<glm::vec3> points2)
	{
		for (unsigned int i = 0; i<points2.size(); i++) {
			points1.push_back(points2[i]);
		}

		return points1;
	}
}