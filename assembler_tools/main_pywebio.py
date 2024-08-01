import os
import tornado.ioloop
import tornado.web
from tornado.web import HTTPError
from typing import Optional

import os
import tornado.ioloop
import tornado.web
import tornado.httpserver
from jinja2 import Environment, FileSystemLoader

# Initialize Jinja2 environment
env = Environment(loader=FileSystemLoader('templates'))


class StaticFileHandler(tornado.web.StaticFileHandler):
    def validate_absolute_path(self, root: str, absolute_path: str) -> Optional[str]:
        root = os.path.abspath(root)
        if not root.endswith(os.path.sep):
            root += os.path.sep
        if not (absolute_path + os.path.sep).startswith(root):
            raise HTTPError(403, "%s is not in root static directory", self.path)
        if os.path.isdir(absolute_path) and self.default_filename is not None:
            if not self.request.path.endswith("/"):
                if self.request.path.startswith("//"):
                    raise HTTPError(
                        403, "cannot redirect path with two initial slashes"
                    )
                self.redirect(self.request.path + "/", permanent=True)
                return None
            absolute_path = os.path.join(absolute_path, self.default_filename)
        if not os.path.exists(absolute_path):
            raise HTTPError(404)
        if not (os.path.isfile(absolute_path) or os.path.isdir(absolute_path)):  # this is the only line changed
            raise HTTPError(403, "%s is not a file", self.path)
        return absolute_path

    async def get(self, path: str, include_body: bool = True) -> None:
        self.path = self.parse_url_path(path)
        del path  # make sure we don't refer to path instead of self.path again
        absolute_path = self.get_absolute_path(self.root, self.path)
        # If the path is a directory, render the index.html
        if os.path.isdir(absolute_path):
            samples = os.listdir(absolute_path)
            return 'test' + ' '.join(samples)
        return await super().get(path, include_body)

    def write_error(self, status_code, **kwargs):
        if status_code == 404:
            self.write("404 Not Found")
        else:
            super().write_error(status_code, **kwargs)


class NotFoundHandler(tornado.web.RequestHandler):
    def prepare(self):
        self.set_status(404)
        self.write("404 Not Found")
        self.finish()


def make_app():
    return tornado.web.Application([
        (r"/(.*)", StaticFileHandler,
         {"path": os.path.join(os.path.dirname(os.path.dirname(__file__)),
                               '/home/thomas/PycharmProjects/assembler-tools/data')
          }),
        (r".*", NotFoundHandler),
    ])


if __name__ == "__main__":
    app = make_app()
    app.listen(8000)
    print("Server running at http://localhost:8000")
    tornado.ioloop.IOLoop.current().start()

# import json
# import os
# import shutil
# import logging
# from typing import List, Type
#
# from assembler_tools.ani_dendrogram import ani_clustermap
# from assembler_tools.utils import load_plugins, get_relative_path, AssemblyFailedException
# from assembler_tools.Assembly import Assembly
# from assembler_tools.AssemblyImporter import AssemblyImporter
#
# from jinja2 import Environment, PackageLoader, select_autoescape
#
# env = Environment(
#     loader=PackageLoader('assembler_tools', 'templates'),
#     autoescape=select_autoescape(['html'])
# )
#
# import tornado.ioloop
# import tornado.web
# from pywebio.platform.tornado import webio_handler
#
#
# # from pywebio.input import input, FLOAT
# # from pywebio.output import put_text
#
# class StaticFileHandler(tornado.web.StaticFileHandler):
#     def parse_url_path(self, url_path):
#         print(url_path)
#         if url_path == "":
#             url_path = "index.html"
#         return url_path
#
#
# class MainHandler(tornado.web.RequestHandler):
#     template_overview = env.get_template('overview.html')
#     template_assemblies = env.get_template('assemblies.html')
#
#     samples: List[str]
#     relpath: str
#     samples_dir: str
#
#     def get(self):
#         self.write(self.template_overview.render(
#             samples=self.samples,
#             relpath=self.relpath
#         ))
#
#
# def run_server(samples_dir: str = 'data-pb', plugin_dir: str = 'plugins-pb', port=8010, address='localhost'):
#     # if plugin_dir is None:
#     #     plugin_dir = os.environ.get('PLUGIN_DIR', './plugins')
#     # assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"
#     #
#     # importers = load_plugins(plugin_dir)
#     #
#     # MainHandler.samples = [s for s in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, s))]
#     # MainHandler.relpath = get_relative_path(f'{samples_dir}/overview.html', samples_dir)
#     # MainHandler.samples_dir = samples_dir
#
#     print(os.path.join(os.path.dirname(__file__)))
#     application = tornado.web.Application([
#         # (r"/", MainHandler),
#         # (r"/bmi", webio_handler(bmi)),  # bmi is the same function as above
#         (r"/(.*)", StaticFileHandler,
#          {"path": os.path.join(os.path.dirname(__file__)), "default_filename": "index.html"}),
#     ])
#     application.listen(port=port, address=address)
#     tornado.ioloop.IOLoop.current().start()
#
#
# def curate_assemblies(assemblies: [Assembly]) -> Assembly:
#     print(f"Curating {len(assemblies)} assemblies")
#     return assemblies[0]
#
#
# def load_assemblies(sample: str, sample_dir: str, importers: List[Type[AssemblyImporter]]) -> ([Assembly], [str]):
#     print(f"Processing {sample}")
#
#     assemblies: [Assembly] = []
#     messages: [str] = []
#     for importer_class in importers:
#         try:
#             importer = importer_class(sample_dir)
#             print(f"  -> loading {sample} with {importer.name}")
#             assembly = importer.load_assembly()
#             assemblies.append(assembly)
#             assembly.pprint()
#         except AssemblyFailedException as e:
#             logging.warning(str(e))
#             messages.append(e)
#
#     if not assemblies:
#         msg = f"Failed to load any assemblies for {sample}"
#         logging.warning(msg)
#         messages.insert(0, msg)
#     else:
#         print(f"Loaded {len(assemblies)} assemblies for {sample}")
#     return assemblies, messages
#
#
# def process_samples(importers: [Type[AssemblyImporter]], samples_dir: str, overview_file: str = None):
#     if overview_file is None:
#         overview_file = f'{samples_dir}/overview.html'
#     samples = [s for s in os.listdir(samples_dir) if os.path.isdir(os.path.join(samples_dir, s))]
#     print(f"Processing {len(samples)} samples in {samples_dir}")
#
#     template_overview.stream(
#         samples=samples,
#         relpath=get_relative_path(overview_file, samples_dir)
#     ).dump(overview_file)
#
#     for sample in samples:
#         sample_dir = os.path.join(samples_dir, sample)
#         if not os.path.isdir(sample_dir):
#             logging.info(f"Skipping {sample} as it is not a directory")
#             continue
#         assemblies, messages = load_assemblies(sample, sample_dir, importers)
#         ani_html = ani_clustermap(
#             assemblies=assemblies,
#             output_file=f'{samples_dir}/pyskani_similarity_matrix.tsv'
#         )
#
#         with open(f"{sample_dir}/assemblies.json", 'w') as f:
#             json.dump({assembly.assembler: assembly.to_json() for assembly in assemblies}, f, indent=2)
#
#         template_assemblies.stream(
#             messages=messages,
#             sample=sample,
#             assemblies=assemblies,
#             ani_html=ani_html,
#         ).dump(f"{sample_dir}/assemblies.html")
#
#         def link(src, dst):
#             if os.path.isfile(dst):
#                 os.remove(dst)
#             os.link(src, dst)
#
#         copy_or_link = link  # shutil.copy
#         # Add dotplot.js
#         copy_or_link(
#             os.path.join('/home/thomas/PycharmProjects/dotplot.js/dotplot.js'),
#             f"{samples_dir}/dotplot.js"
#         )
#         # use link instead
#
#         # Add assemblies.css
#         copy_or_link(
#             os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.css'),
#             f"{samples_dir}/assemblies.css"
#         )
#
#         # Add assemblies.js
#         copy_or_link(
#             os.path.join(os.path.dirname(__file__), 'templates', 'assemblies.js'),
#             f"{samples_dir}/assemblies.js"
#         )
#
#         # final_assembly = curate_assemblies(assemblies)
#         # print(f"Final assembly for {sample}: {final_assembly}")
#
#
# def cli(samples_dir: str = './data', plugin_dir: str = None):
#     LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
#     logging.basicConfig(level=LOGLEVEL)
#
#     if plugin_dir is None:
#         plugin_dir = os.environ.get('PLUGIN_DIR', './plugins')
#     assert os.path.isdir(plugin_dir), f"Plugin directory {plugin_dir} does not exist"
#
#     importers = load_plugins(plugin_dir)
#
#     process_samples(importers, samples_dir)
#
#
# def main():
#     LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
#     logging.basicConfig(level=LOGLEVEL)
#
#     from fire import Fire
#
#     Fire(run_server)
#
#
# if __name__ == '__main__':
#     main()
