// .vitepress/config.mts
import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'


export default withMermaid({
  title: 'FLUE Documentation',
  base: '/FLUE_Lat/',

  head: [
    ['link', { rel: 'icon', href: '/images/lattice-063.png' }]
  ],

  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      fortran: 'fortran-free-form',
      f90: 'fortran-free-form',
      f95: 'fortran-free-form',
      f03: 'fortran-free-form',
      f08: 'fortran-free-form',
      f77: 'fortran-fixed-form',
      F90: 'fortran-fixed-form',
    },
  },

  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      { text: 'API', link: '/api/' },
      { text: 'Installation', link: '/Installation/' },
      { text: 'Applications', link: '/Applications/' },
    ],

    sidebar: {
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],

      '/Installation/': [
        {
          text: 'Installation Guide',
          items: [
            { text: 'Fortran', link: '/Installation/fortran' },
            { text: 'Python', link: '/Installation/python' },
          ],
        },
      ],

      '/Applications/': [
        {
          text: 'Applications',
          items: [
            { text: 'CSSM_to_OQCD', link: '/Applications/CSSM_to_OQCD' },
            { text: 'ILDG_to_OQCD', link: '/Applications/ILDG_to_QCD' },
            { text: 'OQCD_to_ILDG', link: '/Applications/OQCD_to_ILDG' },
            { text: 'UNIT_to_OQCD', link: '/Applications/UNIT_to_OQCD' },
            { text: 'SU2_HKLS_to_CSSM', link: '/Applications/SU2_HKLS_to_CSSM' },
            { text: 'SU2_HKLS_to_NRQ2CD', link: '/Applications/SU2_HKLS_to_NRQ2CD' },
            { text: 'SU3_heatbath', link: '/Applications/SU3_heatbath' },
            { text: 'SU2_heatbath', link: '/Applications/SU2_heatbath' },
            { text: 'OQCD_stoutSmear', link: '/Applications/OQCD_stoutSmear' },
            { text: 'magnetic', link: '/Applications/magnetic' },
            { text: 'superMWE', link: '/Applications/superMWE' },
          ],
        },
      ],
    },

    search: {
      provider: 'local',
    },

    socialLinks: [
      {
        icon: 'github',
        link: 'https://github.com/RJaBi/FLUE_Lat',
      },
    ],
  },

  mermaid: {},
})