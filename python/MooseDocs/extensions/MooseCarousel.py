from markdown.blockprocessors import BlockProcessor
from MooseCommonExtension import MooseCommonExtension
import glob
import re
from markdown.util import etree

class MooseCarousel(BlockProcessor, MooseCommonExtension):
    """
    Markdown extension for showing a bootstrap carousel of images.
    Markdown syntax is:
      !slideshow <options>
          images/intro.png caption=Some caption
          images/more*.png
    Where <options> are key=value pairs.
    See http://getbootstrap.com/javascript/#carousel for allowed options.
    Additionally, "caption" can also be used on the slideshow line to
    set a default caption.
    It is assumed image names will have the same filepath as on the webserver.
    """

    RE = re.compile(r'^!\ ?slideshow(.*)')
    SETTINGS_RE = re.compile("([^\s=]+)=(.*?)(?=(?:\s[^\s=]+=|$))")
    # If there are multiple carousels on the same page then
    # they need to have different ids
    MATCHES_FOUND = 0

    def parseOptions(self, settings_line):
      """
      Parses a string of space seperated key=value pairs.
      This supports having values with spaces in them.
      So something like "key0=foo bar key1=value1"
      is supported.
      Input:
        settings_line[str]: Line to parse
      Returns:
        dict of values that were parsed
      """
      matches = self.SETTINGS_RE.findall(settings_line.strip())
      if not matches:
        return {}
      options = {}
      for entry in matches:
        options[entry[0].strip()] = entry[1].strip()
      return options

    def parseFilenames(self, filenames_block):
      """
      Parse a set of lines with filenames in them and an optional caption.
      Filenames can contain wildcards and glob will be used to expand them.
      Expected input is similar to:
          images/1.png caption=My caption
          images/other*.png
      Input:
        filenames_block[str]: String block to parse
      Return:
        list of dicts. Each dict has keys of "path" which is the filename path
          and "caption" which is the associated caption. Caption will be "" if not
          specified.
      """
      lines = filenames_block.split("\n")
      files = []
      for line in lines:
        sline = line.strip()
        idx = sline.find("caption=")
        if idx >=0 :
          caption = sline[idx+8:].strip()
          fname = sline[:idx].strip()
        else:
          caption = ""
          fname = sline
        new_files = glob.glob(fname)
        if not new_files:
          # If one of the paths is broken then
          # we return an empty list to indicate
          # an error state
          return []
        for f in new_files:
          files.append({"path": f, "caption": caption})
      return files

    def test(self, parent, block):
      """
      Test to see if we should process this block of markdown.
      Inherited from BlockProcessor.
      """
      return self.RE.search(block)

    def run(self, parent, blocks):
      """
      Called when it is determined that we can process this block.
      This will convert the markdown into HTML
      """
      sibling = self.lastChild(parent)
      block = blocks.pop(0)
      m = self.RE.search(block)

      if m:
        # Parse out the options on the slideshow line
        options = m.group(1)
        parsed_options = self.parseOptions(options)
        block = block[m.end() + 1:] # removes the slideshow line

      block, theRest = self.detab(block)

      if m:
        files = block
        div = etree.SubElement(parent, "div")
        filenames = self.parseFilenames(files)
        if not filenames:
            return self.createErrorElement(files, "No matching files found")
        self.createCarousel(parsed_options, div, filenames)
        # We processed this whole block so mark it as done
        block = ""
      else:
        div = sibling

      self.parser.parseChunk(div, block)

      if theRest:
        blocks.insert(0, theRest)

    def createCarousel(self, options, top_div, files):
        """
        Creates the actual HTML required for the carousel to work.
        Input:
          options[dict]: Set on the slidshow line of the markdown
          top_div: div element that will be the carousel
          files[list]: List of dicts with filename paths and associated captions
        """
        carousel_options = {
            "interval": options.get("interval", "5000"),
            "pause": options.get("pause", "hover"),
            "wrap": options.get("wrap", "true"),
            "keyboard": options.get("keyboard", "true"),
            }

        cid = "carousel%s" % self.MATCHES_FOUND
        top_div.set("id", cid)
        top_div.set("class", "carousel slide")
        top_div.set("data-ride", "carousel")
        top_div.set("data-interval", carousel_options["interval"])
        top_div.set("data-pause", carousel_options["pause"])
        top_div.set("data-wrap", carousel_options["wrap"])
        top_div.set("data-keyboard", carousel_options["keyboard"])
        ol = etree.SubElement(top_div, 'ol')
        ol.set("class", "carousel-indicators")
        default_caption = options.get("caption", "")
        for i in range(len(files)):
            li = etree.SubElement(ol, 'li')
            li.set("data-target", "#%s" % cid)
            if i == 0:
                li.set("class", "active")
            li.set("data-slide-to", str(i))
        inner_div = etree.SubElement(top_div, 'div')
        inner_div.set("class", "carousel-inner")
        inner_div.set("role", "listbox")
        for i, f in enumerate(files):
            item_div = etree.SubElement(inner_div, "div")
            active = ""
            if i == 0:
                active = "active"
            item_div.set("class", "item %s" % active)
            img = etree.SubElement(item_div, "img")
            img.set("src", f["path"])
            caption = f["caption"]
            if not caption:
              caption = default_caption

            if caption:
              cap_div = etree.SubElement(item_div, "div")
              cap_div.set("class", "carousel-caption")
              cap_div.text = caption
        self.addControl(top_div, cid, "prev", "Previous")
        self.addControl(top_div, cid, "next", "Next")
        self.MATCHES_FOUND += 1
        return top_div

    def addControl(self, parent, cid, direction, text):
        """
        Utility function to add the left/right controls on the carousel.
        Input:
          parent: parent element
          cid: id of the carousel element
          direction: "prev" or "next" to indicate which direction this control operates on.
          text: The alternate text for screen readers.
        """
        c = etree.SubElement(parent, "a")
        alt_dir = "left"
        chev = "glyphicon-chevron-left"
        if direction == "next":
            alt_dir = "right"
            chev = "glyphicon-chevron-right"

        c.set("class", "%s carousel-control" % alt_dir)
        c.set("href", "#%s" % cid)
        c.set("role", "button")
        c.set("data-slide", direction)
        s = etree.SubElement(c, "span")

        s.set("class", "glyphicon %s" % chev)
        s.set("aria-hidden", "true")
        s = etree.SubElement(c, "span")
        s.set("class", "sr-only")
        s.text = text
